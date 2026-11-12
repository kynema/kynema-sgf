#ifndef TURBINE_EXTERNAL_UTILS_H
#define TURBINE_EXTERNAL_UTILS_H

#include <numbers>
#include <string>
#include "src/actuator/turbine/turbine_utils.H"
#include "src/actuator/actuator_ops.H"
#include "src/actuator/actuator_utils.H"
#include "src/actuator/FLLCOp.H"
#include "src/utilities/IOManager.H"
#include "src/actuator/turbine/external/TurbineExternal.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::actuator::external {

// Note: inline is necessary to avoid duplicate symbol errors
inline void swap_epsilon(vs::Vector& eps)
{
    const auto x = eps.x();
    const auto y = eps.y();
    eps.x() = y;
    eps.y() = x;
}

template <typename datatype>
void make_component_views(datatype& data)
{
    auto& grid = data.grid();
    auto& tdata = data.meta();
    const int num_blades = tdata.num_blades;
    const int num_pts_blade = tdata.num_pts_blade;
    const int num_vel_pts_blade = tdata.num_vel_pts_blade;

    for (int ib = 0; ib < num_blades; ++ib) {
        ComponentView cv;

        const auto start = (ib * num_pts_blade) + 1;
        const auto start_vel = (ib * num_vel_pts_blade) + 1;
        // clang-format off
            cv.pos = ::kynema_sgf::utils::slice(
                grid.pos, start, num_pts_blade);
            cv.force = ::kynema_sgf::utils::slice(
                grid.force, start, num_pts_blade);
            cv.epsilon = ::kynema_sgf::utils::slice(
                grid.epsilon, start, num_pts_blade);
            cv.orientation = ::kynema_sgf::utils::slice(
                grid.orientation, start, num_pts_blade);
            cv.chord = ::kynema_sgf::utils::slice(
                tdata.chord, start, num_pts_blade);
            cv.vel_rel = ::kynema_sgf::utils::slice(
                tdata.vel_rel, start, num_pts_blade);
            cv.vel= ::kynema_sgf::utils::slice(
                grid.vel, start_vel, num_vel_pts_blade);
            cv.vel_pos= ::kynema_sgf::utils::slice(
                grid.vel_pos, start_vel, num_vel_pts_blade);
        // clang-format on

        tdata.blades.emplace_back(cv);
    }
    if (tdata.num_pts_tower > 0) {
        const int num_pts_tower = tdata.num_pts_tower;
        const int ntwr_start = (num_blades * num_pts_blade) + 1;
        auto& cv = tdata.tower;

        cv.pos =
            ::kynema_sgf::utils::slice(grid.pos, ntwr_start, num_pts_tower);
        cv.force =
            ::kynema_sgf::utils::slice(grid.force, ntwr_start, num_pts_tower);
        cv.epsilon =
            ::kynema_sgf::utils::slice(grid.epsilon, ntwr_start, num_pts_tower);
        cv.orientation = ::kynema_sgf::utils::slice(
            grid.orientation, ntwr_start, num_pts_tower);
        cv.chord =
            ::kynema_sgf::utils::slice(tdata.chord, ntwr_start, num_pts_tower);
    }
    {
        auto& cv = tdata.hub;
        cv.pos = ::kynema_sgf::utils::slice(grid.pos, 0, 1);
        cv.force = ::kynema_sgf::utils::slice(grid.force, 0, 1);
        cv.epsilon = ::kynema_sgf::utils::slice(grid.epsilon, 0, 1);
        cv.orientation = ::kynema_sgf::utils::slice(grid.orientation, 0, 1);
        cv.chord = ::kynema_sgf::utils::slice(tdata.chord, 0, 1);
    }
}

/** Abort when a point of a blade or of the tower has no width to spread its
 *  force with
 *
 *  The Gaussian of a point with an epsilon of 0 is 0 everywhere on the mesh:
 *  the force of that point would be left out of the momentum source without
 *  any message.
 */
inline void check_epsilon(
    const ComponentView& cv,
    const int npts,
    const std::string& label,
    const std::string& inputs)
{
    for (int i = 0; i < npts; ++i) {
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            if (!(cv.epsilon[i][n] > 0.0_rt)) {
                amrex::Abort(
                    "Actuator " + label +
                    ": epsilon is not larger than 0 at "
                    "point " +
                    std::to_string(i) +
                    ", so its force would "
                    "not reach the mesh. Set one of " +
                    inputs);
            }
        }
    }
}

template <typename datatype>
void init_epsilon(datatype& data)
{
    auto& tdata = data.meta();

    // Order of the epsilon components. The inputs are read in the order of
    // the documentation, (chord, thickness, span). The Gaussian of a point
    // applies component n of its epsilon along row n of its orientation
    // matrix (ActSrcLineOp). The swap below stores (thickness, chord, span):
    // the first input value acts along row 1, the second along row 0 and
    // the third along row 2 of the orientation sent by the solver.
    // OpenFAST: row 0 is normal to the chord, row 1 is along the chord and
    // row 2 is along the span (blade frame of AeroDyn), so each input value
    // acts in its documented direction. The tower points use the same swap.
    // The FLLC does not use these components: it reads the first value of
    // epsilon (fllc_parse) and of epsilon_min (read_inputs), and component
    // 0 of the swapped eps_chord, which is the second input value of
    // epsilon_chord (CommunicateModelDataOp).
    swap_epsilon(tdata.eps_inp);
    swap_epsilon(tdata.eps_min);
    swap_epsilon(tdata.eps_chord);
    swap_epsilon(tdata.eps_tower);

    {
        const auto& cd = tdata.nacelle_cd;
        const auto& area = tdata.nacelle_area;
        const auto eps =
            std::sqrt(2.0_rt / std::numbers::pi_v<amrex::Real> * cd * area);

        auto& nac_eps = data.grid().epsilon[0];
        nac_eps.x() = amrex::max<amrex::Real>(eps, tdata.eps_min.x());
        nac_eps.y() = amrex::max<amrex::Real>(eps, tdata.eps_min.y());
        nac_eps.z() = amrex::max<amrex::Real>(eps, tdata.eps_min.z());
    }

    for (int ib = 0; ib < tdata.num_blades; ++ib) {
        auto& cv = tdata.blades[ib];

        for (int i = 0; i < tdata.num_pts_blade; ++i) {
            const auto eps_crd = tdata.eps_chord * cv.chord[i];

            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                cv.epsilon[i][n] = amrex::max<amrex::Real>(
                    tdata.eps_min[n], tdata.eps_inp[n], eps_crd[n]);
            }
        }
        check_epsilon(
            cv, tdata.num_pts_blade,
            data.info().label + " blade " + std::to_string(ib),
            "epsilon, epsilon_chord or epsilon_min");
    }
    {
        auto& cv = tdata.tower;
        for (int i = 0; i < tdata.num_pts_tower; ++i) {
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                cv.epsilon[i][n] = amrex::max<amrex::Real>(
                    tdata.eps_min[n], tdata.eps_inp[n], tdata.eps_tower[n]);
            }
        }
        check_epsilon(
            cv, tdata.num_pts_tower, data.info().label + " tower",
            "epsilon, epsilon_tower or epsilon_min");

        // epsilon_tower is one of three lower bounds. Tell the user when it
        // is given and is not the value that the tower points get
        bool tower_eps_unused = false;
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            if ((tdata.eps_tower[n] > 0.0_rt) &&
                (tdata.eps_tower[n] <
                 amrex::max<amrex::Real>(tdata.eps_min[n], tdata.eps_inp[n]))) {
                tower_eps_unused = true;
            }
        }
        if (tower_eps_unused && (tdata.num_pts_tower > 0)) {
            amrex::Print()
                << "\nWARNING: Actuator " << data.info().label
                << ": epsilon_tower is smaller than epsilon or epsilon_min. "
                << "The tower points use the largest of the three in each "
                << "direction, so the value given for epsilon_tower is not "
                << "the one used." << std::endl;
        }
    }
    {
        // Patches and MPI ranks that do not touch the bounding box of the
        // turbine receive no force. The box was sized from the rotor
        // diameter alone (utils::read_inputs), before any epsilon was
        // known: add the distance up to which the Gaussian of a point is
        // not zero
        amrex::Real eps_max = 0.0_rt;
        for (const auto& eps : data.grid().epsilon) {
            eps_max = amrex::max<amrex::Real>(
                eps_max, eps.x(), eps.y(), eps.z());
        }
        const auto reach =
            ::kynema_sgf::actuator::utils::gaussian3d_cutoff() * eps_max;
        auto& bbox = data.info().bound_box;
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            bbox.setLo(n, bbox.lo(n) - reach);
            bbox.setHi(n, bbox.hi(n) + reach);
        }
    }
}

/** Tell the user, once, when the Gaussian of a point reaches beyond the
 *  bounding box of the turbine
 *
 *  The box is made when the turbine is created, from base_position,
 *  hub_height, rotor_diameter and the largest epsilon, and it does not follow
 *  the turbine. Patches outside of it receive no force. The Gaussian of a
 *  point leaves the box when the rotor of the external solver is not where
 *  these inputs put it, or when the turbine moves by more than a quarter of
 *  the rotor radius (a floating platform).
 *
 *  \return Number of points with a Gaussian that leaves the box
 */
template <typename datatype>
int check_bound_box(datatype& data)
{
    const auto& grid = data.grid();
    const auto& bbox = data.info().bound_box;
    const auto cutoff = ::kynema_sgf::actuator::utils::gaussian3d_cutoff();
    // Tolerance for a point that is on a face of the box of the inputs
    const auto tol = 1.0e-8_rt * (bbox.hi(0) - bbox.lo(0));
    const int npts = static_cast<int>(grid.pos.size());
    int num_out = 0;
    amrex::Real dist_max = 0.0_rt;
    for (int i = 0; i < npts; ++i) {
        const auto& eps = grid.epsilon[i];
        const auto reach =
            cutoff * amrex::max<amrex::Real>(eps.x(), eps.y(), eps.z());
        amrex::Real dist = 0.0_rt;
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            dist = amrex::max<amrex::Real>(
                dist, bbox.lo(n) - (grid.pos[i][n] - reach),
                (grid.pos[i][n] + reach) - bbox.hi(n));
        }
        if (dist > tol) {
            ++num_out;
            dist_max = amrex::max<amrex::Real>(dist_max, dist);
        }
    }
    if ((num_out > 0) && !data.meta().bound_box_warned) {
        data.meta().bound_box_warned = true;
        amrex::OutStream()
            << "\nWARNING: Actuator " << data.info().label << ": the force of "
            << num_out << " of " << npts << " points reaches up to "
            << dist_max << " beyond the region reserved for this turbine on "
            << "the mesh, and patches outside of this region receive no "
            << "force. The region is the box that base_position, "
            << "hub_height and rotor_diameter give (1.25 rotor radii around "
            << "the turbine), enlarged by " << cutoff
            << " times the largest epsilon. Check these inputs against the "
            << "model of the external solver; for a turbine that moves this "
            << "far, give a larger rotor_diameter." << std::endl;
    }
    return num_out;
}

template <typename datatype>
void compute_nacelle_force(datatype& data)
{
    if (!data.info().is_root_proc) {
        return;
    }

    const auto& cd = data.meta().nacelle_cd;
    const auto& area = data.meta().nacelle_area;
    const auto& cd_area = cd * area;
    const auto& ext_tdata = data.meta().ext_data;
    const auto& rho = data.meta().density;

    const auto& eps = data.grid().epsilon[0].x();
    // This assumes a static nacelle
    vs::Vector vel{
        ext_tdata.fluid_velocity(0)[0], ext_tdata.fluid_velocity(1)[0],
        ext_tdata.fluid_velocity(2)[0]};
    amrex::Real correction = 0.0_rt;
    if (eps > 0.0_rt) {
        amrex::Real fac =
            1.0_rt -
            (cd_area / (2.0_rt * ::kynema_sgf::utils::two_pi() * eps * eps));
        correction = 1.0_rt / fac;
    }
    amrex::Real coeff =
        0.5_rt * rho * cd_area * vs::mag(vel) * correction * correction;

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        ext_tdata.force(dir)[0] = static_cast<float>(coeff * vel[dir]);
    }
}

template <typename datatype>
void ext_loads(datatype& data)
{
    if (!data.info().is_root_proc) {
        return;
    }

    auto& meta = data.meta();
    auto& tf = data.meta().ext_data;
    meta.ext_ptr->advance_aero_loads(tf.tid_local);
    if (tf.is_solution0) {
        meta.ext_ptr->init_solution(tf.tid_local);
    }

    compute_nacelle_force<datatype>(data);
}

template <typename datatype>
void ext_step(datatype& data)
{
    if (!data.info().is_root_proc) {
        return;
    }

    auto& meta = data.meta();
    auto& tf = data.meta().ext_data;
    meta.ext_ptr->advance_turbine(tf.tid_local);
}

template <typename datatype>
void scatter_data(datatype& data)
{
    if (!data.info().actuator_in_proc) {
        return;
    }

    // Create an MPI transfer buffer that packs all data in one contiguous
    // array. 3 floats for the position vector, 3 floats for the force
    // vector, and 9 floats for the orientation matrix = 15 floats per
    // actuator node.
    const auto dsize = data.grid().pos.size() * 15;
    amrex::Vector<float> buf(dsize);

    // Copy data into MPI send/recv buffer from the OpenFAST data structure.
    // Note, other procs do not have a valid data in those pointers.
    if (data.info().is_root_proc) {
        BL_PROFILE(
            "kynema-sgf::actuator::external::advance_external_turbine::"
            "scatter1");
        const auto& ext_tdata = data.meta().ext_data;
        // Update hub data prior to copy
        data.meta().ext_ptr->get_hub_stats(ext_tdata.tid_local);
        auto it = buf.begin();
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            std::copy(
                ext_tdata.force(dir),
                ext_tdata.force(dir) + ext_tdata.length_force(dir), it);
            std::advance(it, ext_tdata.length_force(dir));
        }
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            std::copy(
                ext_tdata.position_at_force(dir),
                ext_tdata.position_at_force(dir) +
                    ext_tdata.length_position_at_force(dir),
                it);
            std::advance(it, ext_tdata.length_position_at_force(dir));
        }

        // clang-format off
            std::copy(ext_tdata.orientation(),
                      ext_tdata.orientation() + ext_tdata.length_orientation(), it);
        // clang-format on
    }

    // Broadcast data to all influenced procs from the root process
    const auto& procs = data.info().procs;
    const int tag = 1001;
    if (data.info().is_root_proc) {
        BL_PROFILE(
            "kynema-sgf::actuator::external::advance_external_turbine::"
            "scatter2");
        for (const int ip : procs) {
            if (ip == data.info().root_proc) {
                continue;
            }

            amrex::ParallelDescriptor::Send(
                buf.data(), dsize, ip, tag, data.meta().tcomm);
        }
    } else {
        BL_PROFILE(
            "kynema-sgf::actuator::external::advance_external_turbine::"
            "scatter2");
        amrex::ParallelDescriptor::Recv(
            buf.data(), dsize, data.info().root_proc, tag, data.meta().tcomm);
    }

    // Populate the actuator grid data structures with data from the MPI
    // send/recv buffer.
    {
        BL_PROFILE(
            "kynema-sgf::actuator::external::advance_external_turbine::"
            "scatter3");
        const auto& bp = data.info().base_pos;
        auto& grid = data.grid();
        const auto& npts = grid.pos.size();
        const auto& rho = data.meta().density;
        const size_t ifx = 0;
        const size_t ify = ifx + npts;
        const size_t ifz = ify + npts;
        const size_t ipx = ifz + npts;
        const size_t ipy = ipx + npts;
        const size_t ipz = ipy + npts;
        const size_t iori = ipz + npts;

        for (int i = 0; i < npts; ++i) {
            // Aerodynamic force vectors. Flip sign to get force on fluid.
            // Divide by density as the source term computation will
            // multiply by density before adding to momentum equation.
            //
            grid.force[i].x() = -static_cast<amrex::Real>(buf[ifx + i]) / rho;
            grid.force[i].y() = -static_cast<amrex::Real>(buf[ify + i]) / rho;
            grid.force[i].z() = -static_cast<amrex::Real>(buf[ifz + i]) / rho;

            // Position vectors of the actuator nodes. Add shift to base
            // locations.
            grid.pos[i].x() = static_cast<amrex::Real>(buf[ipx + i]) + bp.x();
            grid.pos[i].y() = static_cast<amrex::Real>(buf[ipy + i]) + bp.y();
            grid.pos[i].z() = static_cast<amrex::Real>(buf[ipz + i]) + bp.z();

            // Copy over the orientation matrix
            //
            // Note that we transpose the orientation matrix when copying
            // from external to Kynema-SGF Tensor data structure. This is done
            // so that post-multiplication of vector transforms from global
            // to local reference frame.
            const auto off =
                static_cast<int>(iori) + (i * AMREX_SPACEDIM * AMREX_SPACEDIM);
            for (int j = 0; j < AMREX_SPACEDIM; ++j) {
                for (int k = 0; k < AMREX_SPACEDIM; ++k) {
                    grid.orientation[i][(j * AMREX_SPACEDIM) + k] =
                        static_cast<amrex::Real>(
                            buf[off + j + (k * AMREX_SPACEDIM)]);
                }
            }
        }

        // Extract the rotor center of rotation
        auto& meta = data.meta();
        meta.rot_center = grid.pos[0];

        // Rotor non-rotating reference frame
        const auto xvec = grid.orientation[0].x().unit();
        const auto yvec = vs::Vector::khat() ^ xvec;
        const auto zvec = xvec ^ yvec;
        meta.rotor_frame.rows(xvec, yvec.unit(), zvec.unit());
    }

    if (data.info().is_root_proc) {
        check_bound_box<datatype>(data);
    }
}

} // namespace kynema_sgf::actuator::external

#endif /* TURBINE_EXTERNAL_UTILS_H */
