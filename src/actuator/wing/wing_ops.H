#ifndef WING_OPS_H
#define WING_OPS_H

#include <fstream>
#include "src/actuator/wing/ActuatorWing.H"
#include "src/actuator/actuator_types.H"
#include "src/actuator/actuator_ops.H"
#include "src/actuator/actuator_restart.H"
#include "src/actuator/aero/AirfoilTable.H"
#include "src/actuator/ActParser.H"
#include "src/actuator/FLLCOp.H"
#include "src/utilities/linear_interpolation.H"
#include "src/utilities/io_utils.H"
#include "src/utilities/IOManager.H"
#include "src/core/vs/vector_space.H"
#include "AMReX_REAL.H"

#include <type_traits>

using namespace amrex::literals;

namespace kynema_sgf::actuator {
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
    view.pos = ::kynema_sgf::utils::slice(grid.pos, 0, meta.num_pts);
    view.vel_pos = ::kynema_sgf::utils::slice(grid.vel_pos, 0, meta.num_pts);
    view.force = ::kynema_sgf::utils::slice(grid.force, 0, meta.num_pts);
    view.epsilon = ::kynema_sgf::utils::slice(grid.epsilon, 0, meta.num_pts);
    view.orientation =
        ::kynema_sgf::utils::slice(grid.orientation, 0, meta.num_pts);
    view.chord = ::kynema_sgf::utils::slice(meta.chord, 0, meta.num_pts);
    view.vel_rel = ::kynema_sgf::utils::slice(meta.vel_rel, 0, meta.num_pts);
    view.vel = ::kynema_sgf::utils::slice(grid.vel, 0, meta.num_pts);

    return view;
}

/** Write the FLLC state of a wing into a checkpoint
 *
 *  FLLC state, plus relative velocity and force at the actuator points
 *  (inputs to the next lift). Wing root only; nothing before the FLLC starts.
 */
template <typename T>
void write_checkpoint(typename T::DataType& data, const std::string& chkname)
{
    auto& wdata = data.meta();
    if (!data.info().is_root_proc || !wdata.fllc || !wdata.fllc->initialized) {
        return;
    }

    auto os = restart::create_file(
        chkname, data.info().label, data.sim().time().new_time());
    fllc_write_state(os, *(wdata.fllc));
    restart::write_list(os, "relative_velocity", wdata.component_view.vel_rel);
    restart::write_list(os, "force", wdata.component_view.force);
}

/** Read the FLLC state of a wing from the restart checkpoint
 *
 *  Read on every process that computes the wing forces.
 *
 *  \return False, with nothing changed, when there is no state for this wing
 *  and time, or it is for other points, positions or FLLC inputs
 */
template <typename T>
bool read_checkpoint(typename T::DataType& data, const std::string& chkname)
{
    auto& wdata = data.meta();
    std::ifstream is;
    if (!restart::open_file(
            is, chkname, data.info().label, data.sim().time().current_time())) {
        return false;
    }

    // Read into copies; nothing changes unless the state fits
    const auto npts = wdata.fllc->correction_velocity.size();
    FLLCData fllc(*(wdata.fllc));
    VecList vel_rel(npts);
    VecList force(npts);
    if (!fllc_read_state(is, fllc) ||
        !restart::read_list(is, "relative_velocity", vel_rel) ||
        !restart::read_list(is, "force", force)) {
        return false;
    }

    *(wdata.fllc) = fllc;
    auto& view = wdata.component_view;
    for (int i = 0; i < static_cast<int>(npts); ++i) {
        view.vel_rel[i] = vel_rel[i];
        view.force[i] = force[i];
    }
    return true;
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
    static constexpr bool update_fields = true;
    static constexpr bool prescribe_motion = true;
    static constexpr bool compute_loads = true;
    static constexpr bool advance_model = true;
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
struct UpdateFieldsOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if_t<
        std::is_base_of_v<WingType, ActTrait> &&
        UseDefaultOp<ActTrait, ActSrcLine>::update_fields>>
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
struct PrescribeModelMotionOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if_t<
        std::is_base_of_v<WingType, ActTrait> &&
        UseDefaultOp<ActTrait, ActSrcLine>::prescribe_motion>>
{
    void operator()(typename ActTrait::DataType& data)
    {
        auto& grid = data.grid();
        auto& wdata = data.meta();
        const int npts = wdata.num_pts;
        const auto& time = data.sim().time();

        // Move force points to location at n+1
        // Also, get velocity of actuator from n to n+1
        // (will be unchanged for "none" and "linear" motion)
        actuator::wing::new_wing_position_velocity(
            grid.pos, wdata.vel_tr, npts, time.current_time(), time.new_time(),
            wdata.motion_type, wdata.s_period, wdata.s_vector);

        // Build the local reference frame
        const vs::Vector wspan = wdata.end - wdata.start;

        // Calculate the aoa based on the pitch actuation table if supplied
        // (interpolate pitch angle at step n+1/2) or the initial pitch angle
        if (!wdata.pitch_timetable_file.empty()) {
            const auto time_interp =
                (time.current_time() + time.new_time()) / 2.0_rt;
            wdata.pitch = kynema_sgf::interp::linear_angle(
                wdata.time_table, wdata.pitch_table, time_interp, 360.0_rt);
            amrex::Print() << "Wing Ops: Current wing pitch (at time "
                           << time_interp << " s) is " << wdata.pitch
                           << " degrees.\n";
            // Recalculate transformation matrix and update on grid struct
            auto tmat = vs::quaternion(wspan, wdata.pitch);
            grid.orientation.assign(npts, tmat);
        }
    }
};

template <typename ActTrait>
struct ComputeLoadsOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if_t<
        std::is_base_of_v<WingType, ActTrait> &&
        UseDefaultOp<ActTrait, ActSrcLine>::compute_loads>>
{
    void operator()(typename ActTrait::DataType& data)
    {
        auto& grid = data.grid();
        auto& wdata = data.meta();
        const int npts = wdata.num_pts;
        const auto& dx = wdata.dx;
        const auto& chord = wdata.chord;
        const auto& aflookup = airfoil_lookup<ActTrait>(data);

        // Build the local reference frame
        const vs::Vector wspan = wdata.end - wdata.start;

        // Use the global coord to orient the blade
        // The default is for inflow in the x direction
        const auto blade_y = wspan.unit();
        const auto blade_z = (wdata.blade_x.unit() ^ blade_y).unit();
        // Ensure orthogonality for new reference frame
        const auto blade_x = (blade_y ^ blade_z).unit();

        // Get flags for each coordinate of force and put in vs::Vector
        const vs::Vector force_coord_flags_vec = {
            static_cast<amrex::Real>(wdata.force_coord_flags[0]),
            static_cast<amrex::Real>(wdata.force_coord_flags[1]),
            static_cast<amrex::Real>(wdata.force_coord_flags[2])};

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
                             kynema_sgf::utils::radians(wdata.pitch);

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
                qfac /= static_cast<amrex::Real>(npts);
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
            wdata.aoa[ip] = kynema_sgf::utils::degrees(aoa);
            wdata.cl[ip] = cl;
            wdata.cd[ip] = cd;

            total_lift += lift;
            total_drag += drag;
        }

        wdata.lift = total_lift;
        wdata.drag = total_drag;

        const auto& time = data.sim().time();
        if (wdata.fllc) {
            // Restart: continue from the checkpoint state (else the
            // correction restarts from zero and takes many steps to
            // recover). Without a state, skip the first force call (at the
            // checkpoint time, absent in an unrestarted run) and start at
            // the first step's call.
            auto& io = data.sim().io_manager();
            const bool restart_call =
                wdata.at_first_force_call && io.is_restart();
            bool has_state = false;
            if (restart_call) {
                has_state = restart::has_file(
                    io.restart_file(), data.info().label, time.current_time());
            }

            const bool due =
                !(wdata.fllc->initialized) &&
                (time.current_time() > wdata.fllc->fllc_start_time);
            const bool waits = due && restart_call && !has_state;
            const bool started = due && !waits;
            if (started) {
                wdata.component_view =
                    kynema_sgf::actuator::wing::make_component_view<ActTrait>(
                        data);
                fllc_init(
                    *(wdata.fllc.get()), wdata.component_view,
                    wdata.epsilon_chord[0], data.info().label,
                    data.info().is_root_proc);
            }

            if (started && restart_call) {
                const bool found =
                    kynema_sgf::actuator::wing::read_checkpoint<ActTrait>(
                        data, io.restart_file());
                if (!found && data.info().is_root_proc) {
                    amrex::OutStream()
                        << "\nWARNING: Actuator " << data.info().label
                        << ": the state of the filtered lifting line "
                        << "correction in the checkpoint " << io.restart_file()
                        << " does not fit this wing: it is for another "
                        << "number of points, for other positions of the "
                        << "points or for other inputs of the correction. "
                        << "The correction starts from zero." << std::endl;
                }
            }
            if (waits && data.info().is_root_proc) {
                amrex::OutStream()
                    << "\nWARNING: Actuator " << data.info().label
                    << ": the checkpoint " << io.restart_file()
                    << " has no state of the filtered lifting line "
                    << "correction for this time. If the correction had not "
                    << "started when the checkpoint was written, nothing is "
                    << "lost. If it had (a checkpoint written without the "
                    << "state), the restart loses it. The correction starts "
                    << "from zero in the first time step." << std::endl;
            }
        }
        wdata.at_first_force_call = false;
    }
};

template <typename ActTrait>
struct AdvanceModelOp<
    ActTrait,
    ActSrcLine,
    typename std::enable_if_t<
        std::is_base_of_v<WingType, ActTrait> &&
        UseDefaultOp<ActTrait, ActSrcLine>::advance_model>>
{
    void operator()(typename ActTrait::DataType& /*data*/) {}
};

template <typename ActTrait, typename SrcTrait>
struct ProcessOutputsOp<
    ActTrait,
    SrcTrait,
    typename std::enable_if_t<
        std::is_base_of_v<WingType, ActTrait> &&
        UseDefaultOp<ActTrait, SrcTrait>::process_outputs>>
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
        ioutils::query_with_legacy_name(
            pp, "output_interval", "output_frequency", m_out_freq);
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
} // namespace kynema_sgf::actuator

#endif /* WING_OPS_H */
