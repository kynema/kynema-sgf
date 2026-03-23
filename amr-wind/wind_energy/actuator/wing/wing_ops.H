#ifndef WING_OPS_H
#define WING_OPS_H

#include "amr-wind/wind_energy/actuator/wing/ActuatorWing.H"
#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/wind_energy/actuator/aero/AirfoilTable.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "amr-wind/wind_energy/actuator/FLLCOp.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "amr-wind/core/vs/vector_space.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator {
namespace wing {

// Initialize core data structures when modeling fixed wings
void init_data_structures(WingBaseData& /*wdata*/, ActGrid& /*grid*/);

void prepare_netcdf_file(
    const std::string& /*ncfile*/,
    const WingBaseData& /*meta*/,
    const ActInfo& /*info*/,
    const ActGrid& /*grid*/);

void write_netcdf(
    const std::string& /*ncfile*/,
    const WingBaseData& /*meta*/,
    const ActInfo& /*info*/,
    const ActGrid& /*grid*/,
    amrex::Real /*time*/);

void refresh_wing_position(VecList& vpoints, VecList fpoints, int npts);

void new_wing_position_velocity(
    VecList& points,
    vs::Vector& vtr,
    int npts,
    amrex::Real tn,
    amrex::Real tnp1,
    const std::string& motion,
    amrex::Real period,
    vs::Vector svec);

template <typename T>
ComponentView make_component_view(typename T::DataType& data)
{
    auto& grid = data.grid();
    auto& meta = data.meta();
    ComponentView view;
    view.pos = ::amr_wind::utils::slice(grid.pos, 0, meta.num_pts);
    view.vel_pos = ::amr_wind::utils::slice(grid.vel_pos, 0, meta.num_pts);
    view.force = ::amr_wind::utils::slice(grid.force, 0, meta.num_pts);
    view.epsilon = ::amr_wind::utils::slice(grid.epsilon, 0, meta.num_pts);
    view.orientation =
        ::amr_wind::utils::slice(grid.orientation, 0, meta.num_pts);
    view.chord = ::amr_wind::utils::slice(meta.chord, 0, meta.num_pts);
    view.vel_rel = ::amr_wind::utils::slice(meta.vel_rel, 0, meta.num_pts);
    view.vel = ::amr_wind::utils::slice(grid.vel, 0, meta.num_pts);

    return view;
}

} // namespace wing

namespace ops {

template <typename T>
const typename AirfoilTraits<T>::AirfoilLookup&
airfoil_lookup(typename T::DataType& data);

template <typename ActTrait, typename SrcTrait>
struct UseDefaultOp
{
    static constexpr bool update_pos = true;
    static constexpr bool update_vel = true;
    static constexpr bool compute_force = true;
    static constexpr bool process_outputs = true;
};

template <typename ActTrait>
struct UpdatePosOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if_t<
        std::is_base_of_v<WingType, ActTrait> &&
        UseDefaultOp<ActTrait, ActSrcLine>::update_pos>>
{
    void operator()(typename ActTrait::DataType& data)
    {
        // Put wing at current (time n) position from force points
        auto& grid = data.grid();
        const int npts = data.meta().num_pts;
        actuator::wing::refresh_wing_position(grid.vel_pos, grid.pos, npts);
    }
};

template <typename ActTrait>
struct UpdateVelOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if_t<
        std::is_base_of_v<WingType, ActTrait> &&
        UseDefaultOp<ActTrait, ActSrcLine>::update_vel>>
{

    void operator()(typename ActTrait::DataType& data)
    {
        auto& meta = data.meta();
        if (meta.fllc) {
            FLLCOp()(meta.component_view, *(meta.fllc));
        }
    }
};

template <typename ActTrait>
struct ComputeForceOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if_t<
        std::is_base_of_v<WingType, ActTrait> &&
        UseDefaultOp<ActTrait, ActSrcLine>::compute_force>>
{
    void operator()(typename ActTrait::DataType& data)
    {
        auto& grid = data.grid();
        auto& wdata = data.meta();
        const int npts = wdata.num_pts;
        const auto& dx = wdata.dx;
        const auto& chord = wdata.chord;
        const auto& aflookup = airfoil_lookup<ActTrait>(data);
        const auto& time = data.sim().time();

        // Move force points to location at n+1
        // Also, get velocity of actuator from n to n+1
        // (will be unchanged for "none" and "linear" motion)
        actuator::wing::new_wing_position_velocity(
            grid.pos, wdata.vel_tr, npts, time.current_time(), time.new_time(),
            wdata.motion_type, wdata.s_period, wdata.s_vector);

        // Build the local reference frame
        const vs::Vector wspan = wdata.end - wdata.start;

        // Use the global coord to orient the blade
        // The default is for inflow in the x direction
        const auto blade_y = wspan.unit();
        const auto blade_z = (wdata.blade_x.unit() ^ blade_y).unit();
        // Ensure orthogonality for new reference frame
        const auto blade_x = (blade_y ^ blade_z).unit();

        // Calculate the aoa based on the pitch actuation table if supplied
        // (interpolate pitch angle at step n+1/2) or the initial pitch angle
        if (!wdata.pitch_timetable_file.empty()) {
            const auto time_interp =
                (time.current_time() + time.new_time()) / 2.0_rt;
            wdata.pitch = amr_wind::interp::linear_angle(
                wdata.time_table, wdata.pitch_table, time_interp, 360.0_rt);
            amrex::Print() << "Wing Ops: Current wing pitch (at time "
                           << time_interp << " s) is " << wdata.pitch
                           << " degrees.\n";
            // Recalculate transformation matrix and update on grid struct
            auto tmat = vs::quaternion(wspan, wdata.pitch);
            grid.orientation.assign(npts, tmat);
        }

        // Get flags for each coordinate of force and put in vs::Vector
        const vs::Vector force_coord_flags_vec = {
            (amrex::Real)wdata.force_coord_flags[0],
            (amrex::Real)wdata.force_coord_flags[1],
            (amrex::Real)wdata.force_coord_flags[2]};

        // Calculate the local force using sampled velocity (at n)
        amrex::Real total_lift = 0.0_rt;
        amrex::Real total_drag = 0.0_rt;
        for (int ip = 0; ip < npts; ++ip) {

            // Wind vector is relative to actuator motion
            vs::Vector windvector;
            windvector[0] = (grid.vel[ip] - wdata.vel_tr) & blade_x;
            windvector[1] = 0;
            windvector[2] = (grid.vel[ip] - wdata.vel_tr) & blade_z;

            const auto vmag = vs::mag(windvector);
            const auto aoa = std::atan2(windvector[2], windvector[0]) +
                             amr_wind::utils::radians(wdata.pitch);

            // Get Cl, Cd values
            amrex::Real cl, cd;
            aflookup(aoa, cl, cd);

            // Calculate factor for qval: 0.5_rt * Uinf * Uinf * dx
            // but replace velocity magnitude if specified
            auto qfac =
                (wdata.prescribed_uinf < 0.0_rt)
                    ? 0.5_rt * vmag * vmag
                    : 0.5_rt * wdata.prescribed_uinf * wdata.prescribed_uinf;
            // Length scale in spanwise direction
            if (wdata.gauss_2D && wdata.normalize_2D_spanwise) {
                // Do not multiply by a length scale, but divide by number of
                // points because the Gaussians of each point all overlap
                qfac /= (amrex::Real)npts;
                // Convert 3D Gaussian factor to 2D factor
                constexpr amrex::Real G2D_from_3D =
                    std::numbers::inv_pi_v<amrex::Real> /
                    0.17958712212516656_rt;
                qfac *= G2D_from_3D;
                // Get rid of epsilon in spanwise direction
                const amrex::Real eps_span = grid.epsilon[ip][1];
                qfac *= eps_span;
            } else {
                // 3D implementation: multiply by spanwise dx
                qfac *= dx[ip];
            }

            const auto qval = qfac * chord[ip];
            const auto lift = qval * cl;
            const auto drag = qval * cd;

            // Determine unit vector parallel and perpendicular to velocity
            const auto drag_dir =
                ((blade_x * windvector.x()) + (blade_z * windvector.z()))
                    .unit();
            const auto lift_dir = (drag_dir ^ blade_y).unit();

            // Compute force on fluid from this section of wing
            grid.force[ip] = -((lift_dir * lift) + (drag * drag_dir));

            // In global coords, zero disabled force components
            grid.force[ip] *= force_coord_flags_vec;

            // Assign values for output
            wdata.vel_rel[ip] = windvector;
            wdata.aoa[ip] = amr_wind::utils::degrees(aoa);
            wdata.cl[ip] = cl;
            wdata.cd[ip] = cd;

            total_lift += lift;
            total_drag += drag;
        }

        wdata.lift = total_lift;
        wdata.drag = total_drag;

        if (wdata.fllc) {
            if (!(wdata.fllc->initialized) &&
                (time.current_time() > wdata.fllc->fllc_start_time)) {
                wdata.component_view =
                    amr_wind::actuator::wing::make_component_view<ActTrait>(
                        data);
                fllc_init(
                    *(wdata.fllc.get()), wdata.component_view,
                    wdata.epsilon_chord[0]);
            }
        }
    }
};

template <typename ActTrait, typename SrcTrait>
struct ProcessOutputsOp<
    ActTrait,
    SrcTrait,
    typename std::enable_if_t<
        std::is_base_of_v<WingType, ActTrait> &&
        UseDefaultOp<ActTrait, ActSrcLine>::process_outputs>>
{
private:
    typename ActTrait::DataType& m_data;

    std::string m_out_dir;
    std::string m_nc_filename;
    int m_out_freq{10};

public:
    explicit ProcessOutputsOp(typename ActTrait::DataType& data) : m_data(data)
    {}

    void read_io_options(const utils::ActParser& pp)
    {
        pp.query("output_frequency", m_out_freq);
    }

    void prepare_outputs(const std::string& out_dir)
    {
        m_nc_filename = out_dir + "/" + m_data.info().label + ".nc";
        wing::prepare_netcdf_file(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid());
    }

    void write_outputs()
    {
        const auto& time = m_data.sim().time();
        const int tidx = time.time_index();
        if (tidx % m_out_freq != 0) {
            return;
        }

        wing::write_netcdf(
            m_nc_filename, m_data.meta(), m_data.info(), m_data.grid(),
            time.new_time());
    }
};

} // namespace ops
} // namespace amr_wind::actuator

#endif /* WING_OPS_H */
