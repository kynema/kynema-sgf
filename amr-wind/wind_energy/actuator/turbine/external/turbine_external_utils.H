#ifndef TURBINE_EXTERNAL_UTILS_H
#define TURBINE_EXTERNAL_UTILS_H

#include <numbers>
#include "amr-wind/wind_energy/actuator/turbine/turbine_utils.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/wind_energy/actuator/actuator_utils.H"
#include "amr-wind/wind_energy/actuator/FLLCOp.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/wind_energy/actuator/turbine/external/TurbineExternal.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator::external {

void swap_epsilon(vs::Vector& eps)
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

        // This might be another problem! for Kynema vs OpenFAST
        // Only relevant for disks, I think. Would need to adjust start_vel
        const auto start = ib * num_pts_blade + 1;
        const auto start_vel = ib * num_vel_pts_blade;
        // clang-format off
            cv.pos = ::amr_wind::utils::slice(
                grid.pos, start, num_pts_blade);
            cv.force = ::amr_wind::utils::slice(
                grid.force, start, num_pts_blade);
            cv.epsilon = ::amr_wind::utils::slice(
                grid.epsilon, start, num_pts_blade);
            cv.orientation = ::amr_wind::utils::slice(
                grid.orientation, start, num_pts_blade);
            cv.chord = ::amr_wind::utils::slice(
                tdata.chord, start, num_pts_blade);
            cv.vel_rel = ::amr_wind::utils::slice(
                tdata.vel_rel, start, num_pts_blade);
            cv.vel= ::amr_wind::utils::slice(
                grid.vel, start_vel, num_vel_pts_blade);
            cv.vel_pos= ::amr_wind::utils::slice(
                grid.vel_pos, start_vel, num_vel_pts_blade);
        // clang-format on

        tdata.blades.emplace_back(cv);
    }
    if (tdata.num_pts_tower > 0) {
        const int num_pts_tower = tdata.num_pts_tower;
        const int ntwr_start = num_blades * num_pts_blade + 1;
        auto& cv = tdata.tower;

        cv.pos = ::amr_wind::utils::slice(grid.pos, ntwr_start, num_pts_tower);
        cv.force =
            ::amr_wind::utils::slice(grid.force, ntwr_start, num_pts_tower);
        cv.epsilon =
            ::amr_wind::utils::slice(grid.epsilon, ntwr_start, num_pts_tower);
        cv.orientation = ::amr_wind::utils::slice(
            grid.orientation, ntwr_start, num_pts_tower);
        cv.chord =
            ::amr_wind::utils::slice(tdata.chord, ntwr_start, num_pts_tower);
    }
    {
        auto& cv = tdata.hub;
        cv.pos = ::amr_wind::utils::slice(grid.pos, 0, 1);
        cv.force = ::amr_wind::utils::slice(grid.force, 0, 1);
        cv.epsilon = ::amr_wind::utils::slice(grid.epsilon, 0, 1);
        cv.orientation = ::amr_wind::utils::slice(grid.orientation, 0, 1);
        cv.chord = ::amr_wind::utils::slice(tdata.chord, 0, 1);
    }
}

template <typename datatype>
void init_epsilon(datatype& data)
{
    auto& tdata = data.meta();

    // Swap order of epsilon based on turbine orientation
    // Input order is (chord, span, thickness)
    // Output order should depend on turbine type; currently ...
    // -- this function is correct for Kynema
    // -- this function is erroneous for OpenFAST but has been in use
    // (Does not matter when epsilon is uniform)
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
    }
    {
        auto& cv = tdata.tower;
        for (int i = 0; i < tdata.num_pts_tower; ++i) {
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                cv.epsilon[i][n] = amrex::max<amrex::Real>(
                    tdata.eps_min[n], tdata.eps_inp[n], tdata.eps_tower[n]);
            }
        }
    }
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
            (cd_area) / (2.0_rt * ::amr_wind::utils::two_pi() * eps * eps);
        correction = 1.0_rt / fac;
    }
    amrex::Real coeff =
        0.5_rt * rho * cd_area * vs::mag(vel) * correction * correction;

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        ext_tdata.force(dir)[0] = static_cast<float>(coeff * vel[dir]);
    }
}

template <typename datatype>
void ext_step(datatype& data)
{
    if (!data.info().is_root_proc) {
        return;
    }

    auto& meta = data.meta();
    auto& tf = data.meta().ext_data;
    if (tf.is_solution0) {
        meta.ext_ptr->init_solution(tf.tid_local);
    } else {
        meta.ext_ptr->advance_turbine(tf.tid_local);
    }

    meta.ext_ptr->get_hub_stats(tf.tid_local);

    // Populate nacelle force into the OpenFAST data structure so that it
    // gets broadcasted to all influenced processes in subsequent scattering
    // of data.
    compute_nacelle_force<datatype>(data);
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
        BL_PROFILE("amr-wind::actuator::external::compute_force_op::scatter1");
        const auto& ext_tdata = data.meta().ext_data;
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
        BL_PROFILE("amr-wind::actuator::external::compute_force_op::scatter2");
        for (const int ip : procs) {
            if (ip == data.info().root_proc) {
                continue;
            }

            amrex::ParallelDescriptor::Send(
                buf.data(), dsize, ip, tag, data.meta().tcomm);
        }
    } else {
        BL_PROFILE("amr-wind::actuator::external::compute_force_op::scatter2");
        amrex::ParallelDescriptor::Recv(
            buf.data(), dsize, data.info().root_proc, tag, data.meta().tcomm);
    }

    // Populate the actuator grid data structures with data from the MPI
    // send/recv buffer.
    {
        BL_PROFILE("amr-wind::actuator::external::compute_force_op::scatter3");
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
            // from external to AMR-Wind Tensor data structure. This is done
            // so that post-multiplication of vector transforms from global
            // to local reference frame.
            const auto off =
                static_cast<int>(iori) + i * AMREX_SPACEDIM * AMREX_SPACEDIM;
            for (int j = 0; j < AMREX_SPACEDIM; ++j) {
                for (int k = 0; k < AMREX_SPACEDIM; ++k) {
                    grid.orientation[i][j * AMREX_SPACEDIM + k] =
                        static_cast<amrex::Real>(
                            buf[off + j + k * AMREX_SPACEDIM]);
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
}

} // namespace amr_wind::actuator::external

#endif /* TURBINE_EXTERNAL_UTILS_H */
