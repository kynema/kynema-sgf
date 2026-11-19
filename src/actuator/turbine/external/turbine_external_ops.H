#ifndef TURBINE_EXTERNAL_OPS_H
#define TURBINE_EXTERNAL_OPS_H

#include <fstream>
#include <type_traits>
#include <vector>
#include "src/actuator/turbine/turbine_utils.H"
#include "src/actuator/actuator_ops.H"
#include "src/actuator/actuator_utils.H"
#include "src/actuator/actuator_restart.H"
#include "src/actuator/FLLCOp.H"
#include "src/utilities/IOManager.H"
#include "src/actuator/turbine/external/turbine_external_utils.H"

using namespace kynema_sgf::actuator::external;

namespace kynema_sgf::actuator::ops {

template <typename datatype>
void read_ops(datatype& data, const utils::ActParser& pp)
{
    // Tell reading utility the name of the solver. ext_ptr is not set yet
    // (external_determine_root_proc sets it, on the root process only), so
    // the static function is called through the type
    data.meta().solver_name =
        std::remove_pointer_t<decltype(data.meta().ext_ptr)>::identifier();
    // Data common to any turbine actuator simulation
    utils::read_inputs(data.meta(), data.info(), pp);

    auto& tdata = data.meta();

    // Get density value for normalization
    pp.get("density", tdata.density);
    // The forces of the external solver are divided by it
    if (!(tdata.density > 0.0_rt)) {
        amrex::Abort(
            "Actuator " + data.info().label +
            ": density must be larger than 0. Value read: " +
            std::to_string(tdata.density));
    }

    // Copy data to external data holder
    const auto& tinfo = data.info();
    auto& tf = data.meta().ext_data;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        tf.base_pos[i] = static_cast<float>(tinfo.base_pos[i]);
    }

    tf.tlabel = tinfo.label;
    tf.tid_global = tinfo.id;
    tf.num_blades = tdata.num_blades;
    tf.num_pts_blade = tdata.num_pts_blade;
    tf.num_pts_tower = tdata.num_pts_tower;
    tf.dt_cfd = data.sim().time().delta_t();
}

template <typename datatype>
void external_determine_influenced_procs(datatype& data)
{
    auto& info = data.info();
    info.procs =
        utils::determine_influenced_procs(data.sim().mesh(), info.bound_box);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        info.root_proc > -1,
        "info.root_proc = " + std::to_string(info.root_proc));
    // During regrid, the influenced processes might have changed and might
    // no longer include the root proc. We insert it back to ensure that it
    // is always present on the list.
    info.procs.insert(info.root_proc);

    const int iproc = amrex::ParallelDescriptor::MyProc();
    auto in_proc = info.procs.find(iproc);
    info.actuator_in_proc = (in_proc != info.procs.end());
    info.sample_vel_in_proc = info.is_root_proc;
}

template <typename datatype, typename SolverTurbine, typename SolverData>
void external_determine_root_proc(
    datatype& data, amrex::Vector<int>& act_proc_count)
{
    namespace utils = ::kynema_sgf::actuator::utils;
    auto& info = data.info();
    info.procs =
        utils::determine_influenced_procs(data.sim().mesh(), info.bound_box);

    utils::determine_root_proc(info, act_proc_count);

    // TODO: This function is doing a lot more than advertised by the name.
    // Should figure out a better way to perform the extra work.

    // For OpenFAST we only need velocities sampled in root process
    info.sample_vel_in_proc = info.is_root_proc;

    // Initialize the OpenFAST object and register this turbine in the root
    // process
    if (info.is_root_proc) {
        auto& tdata = data.meta();
        auto& ext_mgr = data.sim().ext_solver_manager();
        ext_mgr.create(::ext_turb::ext_id<SolverData>(), data.sim());
        tdata.ext_ptr =
            &(ext_mgr.template get<
                ::ext_turb::ExtTurbIface<SolverTurbine, SolverData>>());
        tdata.ext_ptr->register_turbine(tdata.ext_data);
    }
}

template <typename datatype>
void init_data_op(datatype& data)
{

    const auto& info = data.info();
    auto& tdata = data.meta();

    // Initialize our communicator for broadcasting data
    amrex::ParallelDescriptor::Comm_dup(
        amrex::ParallelDescriptor::Communicator(), tdata.tcomm);

    amrex::Array<int, 5> sz_info = {0, 0, 0, 0, 0};
    if (info.is_root_proc) {
        tdata.ext_ptr->init_turbine(tdata.ext_data.tid_local);

        const auto& tf = tdata.ext_data;
        sz_info[0] = tf.num_blades;
        sz_info[1] = tf.length_force(0);
        sz_info[2] = tf.length_fluid_velocity(0);
        sz_info[3] = tf.num_pts_tower;
        sz_info[4] = tf.num_vel_pts_blade();
    }

    // Broadcast data to everyone
    amrex::ParallelDescriptor::Bcast(
        sz_info.begin(), sz_info.size(), info.root_proc, tdata.tcomm);

    {
        tdata.num_blades = sz_info[0];
        // FLLC data sized from the num_blades input: one per blade of the
        // external solver
        if (!tdata.fllc.empty()) {
            tdata.fllc.resize(tdata.num_blades, tdata.fllc.front());
        }
        // back calculate what the value per blade for number of points in
        // the openfast data structure
        tdata.num_vel_pts_blade = sz_info[4];
        data.grid().resize(sz_info[1], sz_info[2]);
        tdata.chord.resize(sz_info[1]);
        tdata.num_pts_tower = sz_info[3];
    }

    tdata.vel_rel.assign(sz_info[1], vs::Vector::zero());

    if (info.is_root_proc) {
        // copy chord data
        int npts = sz_info[1];
        for (int i = 0; i < npts; ++i) {
            tdata.chord[i] =
                static_cast<amrex::Real>(tdata.ext_data.chord_at_force()[i]);
        }
    }

    amrex::ParallelDescriptor::Bcast(
        tdata.chord.data(), tdata.chord.size(), info.root_proc, tdata.tcomm);

    make_component_views<datatype>(data);
    init_epsilon<datatype>(data);
}

/** \brief Copy velocity points of the external solver
 *
 *  Copies vel_pos.size() points, from point start of the solver's list of
 *  velocity points on, shifted by the base position of the turbine. The
 *  solver holds the points of its current time, and only on the root
 *  process of the turbine.
 */
template <typename datatype, typename VecType>
void copy_ext_vel_pos(const datatype& data, VecType& vel_pos, const int start)
{
    const auto& tdata = data.meta();
    const auto& bp = data.info().base_pos;
    const auto& pxvel = tdata.ext_data.position_at_vel(0);
    const auto& pyvel = tdata.ext_data.position_at_vel(1);
    const auto& pzvel = tdata.ext_data.position_at_vel(2);
    for (int i = 0; i < static_cast<int>(vel_pos.size()); ++i) {
        vel_pos[i].x() = static_cast<amrex::Real>(pxvel[start + i]) + bp.x();
        vel_pos[i].y() = static_cast<amrex::Real>(pyvel[start + i]) + bp.y();
        vel_pos[i].z() = static_cast<amrex::Real>(pzvel[start + i]) + bp.z();
    }
}

template <typename datatype>
void update_pos_op(datatype& data)
{
    // Return early if this is not the root process for this turbine
    //
    // This is handled in Actuator class, but we add a check here just as a
    // safeguard
    if (!data.info().is_root_proc) {
        return;
    }

    copy_ext_vel_pos(data, data.grid().vel_pos, 0);
}

template <typename datatype>
void update_fields_op(datatype& data)
{
    // Return early if this is not the root process for this turbine
    //
    // This is handled in Actuator class, but we add a check here just as a
    // safeguard
    if (!data.info().is_root_proc) {
        return;
    }
    auto& tdata = data.meta();

    const auto& uvel = tdata.ext_data.fluid_velocity(0);
    const auto& vvel = tdata.ext_data.fluid_velocity(1);
    const auto& wvel = tdata.ext_data.fluid_velocity(2);
    const auto& vel = data.grid().vel;

    if (!tdata.fllc.empty()) {
        namespace interp = ::kynema_sgf::interp;

        // The relative velocity needed for the FLLC is the fluid velocity
        // last sent to the external solver minus the velocity of the blade.
        // The first is known at the velocity points and the second at the
        // force points, each with its own numbering: the hub first, then
        // the points of every blade.
        const auto& xdot = tdata.ext_data.solid_velocity(0);
        const auto& ydot = tdata.ext_data.solid_velocity(1);
        const auto& zdot = tdata.ext_data.solid_velocity(2);
        VecList fluid_vel(tdata.num_vel_pts_blade);
        // Loop through each blade and apply the FLLC
        for (int ib = 0; ib < tdata.num_blades; ++ib) {
            auto& blade = tdata.blades[ib];
            auto& fllc = tdata.fllc[ib];
            if (!fllc.initialized) {
                continue;
            }
            // Same first points as the blade views (make_component_views)
            const int start = (ib * tdata.num_pts_blade) + 1;
            const int start_vel = (ib * tdata.num_vel_pts_blade) + 1;
            for (int i = 0; i < tdata.num_vel_pts_blade; ++i) {
                fluid_vel[i][0] = uvel[start_vel + i];
                fluid_vel[i][1] = vvel[start_vel + i];
                fluid_vel[i][2] = wvel[start_vel + i];
            }
            for (int i = 0; i < tdata.num_pts_blade; ++i) {
                const vs::Vector solid_vel{
                    static_cast<amrex::Real>(xdot[start + i]),
                    static_cast<amrex::Real>(ydot[start + i]),
                    static_cast<amrex::Real>(zdot[start + i])};
                if (fllc.mismatched_points) {
                    blade.vel_rel[i] = interp::linear(
                                           fllc.span_distance_vel, fluid_vel,
                                           fllc.span_distance_force[i]) -
                                       solid_vel;
                } else {
                    blade.vel_rel[i] = fluid_vel[i] - solid_vel;
                }
            }
            FLLCOp()(blade, fllc);
        }
    }

    for (int i = 0; i < tdata.ext_data.length_fluid_velocity(0); ++i) {
        uvel[i] = static_cast<float>(vel[i].x());
        vvel[i] = static_cast<float>(vel[i].y());
        wvel[i] = static_cast<float>(vel[i].z());
    }
}

/** Write the external solver checkpoint with a CFD checkpoint
 *
 *  Called on every rank at each CFD checkpoint, at the end of a time step
 *  (turbine at the new time since pre_advance_work). Turbine root only; no
 *  collectives.
 */
template <typename datatype>
void external_write_solver_checkpoint(datatype& data)
{
    if (!data.info().is_root_proc) {
        return;
    }
    auto& tdata = data.meta();
    tdata.ext_ptr->write_turbine_checkpoint(tdata.ext_data.tid_local);
}

/** Restart data of an external turbine, with each CFD checkpoint
 *  (ops::write_checkpoint of TurbineFast and TurbineKynemaFMB)
 *
 *  External solver checkpoint; once the FLLC has started, per blade its
 *  state and the fluid velocity last sent to the external solver. Turbine
 *  root only.
 */
template <typename datatype>
void external_write_checkpoint(datatype& data, const std::string& chkname)
{
    external_write_solver_checkpoint(data);

    auto& tdata = data.meta();
    if (!data.info().is_root_proc || tdata.fllc.empty() ||
        !tdata.fllc.front().initialized) {
        return;
    }

    auto os = restart::create_file(
        chkname, data.info().label, data.sim().time().new_time());
    restart::write_value(os, "num_blades", tdata.num_blades);

    const auto& uvel = tdata.ext_data.fluid_velocity(0);
    const auto& vvel = tdata.ext_data.fluid_velocity(1);
    const auto& wvel = tdata.ext_data.fluid_velocity(2);
    VecList fluid_vel(tdata.num_vel_pts_blade);
    for (int ib = 0; ib < tdata.num_blades; ++ib) {
        fllc_write_state(os, tdata.fllc[ib]);
        // Same first point as the blade views (make_component_views)
        const int start_vel = (ib * tdata.num_vel_pts_blade) + 1;
        for (int i = 0; i < tdata.num_vel_pts_blade; ++i) {
            fluid_vel[i][0] = uvel[start_vel + i];
            fluid_vel[i][1] = vvel[start_vel + i];
            fluid_vel[i][2] = wvel[start_vel + i];
        }
        restart::write_list(os, "fluid_velocity", fluid_vel);
    }
}

/** Take the state of the filtered lifting line correction of an external
 *  turbine from the checkpoint that the run restarts from
 *
 *  \return False when the checkpoint has no state of this turbine for the
 *  current time, or when the state is for other numbers of blades or points,
 *  for other positions of the points or for other inputs of the correction.
 *  Nothing is changed then.
 */
template <typename datatype>
bool external_read_checkpoint(datatype& data, const std::string& chkname)
{
    auto& tdata = data.meta();
    std::ifstream is;
    if (!restart::open_file(
            is, chkname, data.info().label, data.sim().time().current_time())) {
        return false;
    }
    int num_blades = 0;
    if (!restart::read_value(is, "num_blades", num_blades) ||
        (num_blades != tdata.num_blades)) {
        return false;
    }

    // Read into copies, so that nothing is changed when the state of one
    // of the blades does not fit
    std::vector<FLLCData> fllc(tdata.fllc);
    std::vector<VecList> fluid_vel(num_blades);
    for (int ib = 0; ib < num_blades; ++ib) {
        fluid_vel[ib].resize(tdata.num_vel_pts_blade);
        if (!fllc_read_state(is, fllc[ib]) ||
            !restart::read_list(is, "fluid_velocity", fluid_vel[ib])) {
            return false;
        }
    }

    tdata.fllc = fllc;
    const auto& uvel = tdata.ext_data.fluid_velocity(0);
    const auto& vvel = tdata.ext_data.fluid_velocity(1);
    const auto& wvel = tdata.ext_data.fluid_velocity(2);
    for (int ib = 0; ib < num_blades; ++ib) {
        const int start_vel = (ib * tdata.num_vel_pts_blade) + 1;
        for (int i = 0; i < tdata.num_vel_pts_blade; ++i) {
            uvel[start_vel + i] = static_cast<float>(fluid_vel[ib][i][0]);
            vvel[start_vel + i] = static_cast<float>(fluid_vel[ib][i][1]);
            wvel[start_vel + i] = static_cast<float>(fluid_vel[ib][i][2]);
        }
    }
    return true;
}

template <typename datatype>
void update_external_turbine_loads(datatype& data)
{
    // Update external turbine loads
    ext_loads<datatype>(data);
}

template <typename datatype>
void advance_external_turbine(datatype& data)
{
    // Advance external solver by specified number of sub-steps
    ext_step<datatype>(data);
}

template <typename SolverTurbine, typename SolverData, typename SrcTrait>
struct CommunicateModelDataOp<
    TurbineFromExtSolver<SolverTurbine, SolverData>,
    SrcTrait>
{
    using DataType =
        typename TurbineFromExtSolver<SolverTurbine, SolverData>::DataType;

    void operator()(DataType& data)
    {
        BL_PROFILE(
            "kynema-sgf::actuator::CommunicateModelDataOp<"
            "TurbineFromExtSolver>");
        scatter_data<DataType>(data);

        const auto& time = data.sim().time();

        auto& tdata = data.meta();
        if (!tdata.fllc.empty()) {
            // A run that restarts with the correction on goes on from the state
            // of the checkpoint. Without it the correction starts from zero and
            // needs many time steps to come back. The checkpoint holds the
            // state only when the correction had started before it was written.
            // When it holds none, the correction does not start at the first
            // force call of the restart: that call is made at the time of the
            // checkpoint, where a run without restart makes none. It starts at
            // the force call of the first time step, as it does in that run.
            // The state belongs to the history of the turbine: when the
            // external solver starts the turbine anew (init) on a restart of
            // the flow, the state is not used and the correction starts as
            // without one.
            auto& io = data.sim().io_manager();
            const bool restart_call =
                tdata.at_first_force_call && io.is_restart();
            const bool solver_starts_over =
                tdata.ext_data.sim_mode == ::ext_turb::SimMode::init;
            bool state_in_checkpoint = false;
            if (restart_call) {
                state_in_checkpoint = restart::has_file(
                    io.restart_file(), data.info().label, time.current_time());
            }
            const bool has_state = state_in_checkpoint && !solver_starts_over;

            bool started = false;
            bool waits = false;
            for (int i = 0; i < tdata.num_blades; ++i) {
                if (!(tdata.fllc[i].initialized) &&
                    (time.current_time() > tdata.fllc[i].fllc_start_time)) {
                    if (restart_call && !has_state) {
                        waits = true;
                        continue;
                    }
                    // fllc_init keeps point distances from the first point,
                    // so both point sets must be of one instant. After an
                    // external step the grid has new force points but old
                    // velocity points (update_pos_op ran before the step):
                    // take the new velocity points from the solver. Root only.
                    ComponentView view = tdata.blades[i];
                    VecList vel_pos_now(tdata.num_vel_pts_blade);
                    if (data.info().is_root_proc) {
                        copy_ext_vel_pos(
                            data, vel_pos_now,
                            (i * tdata.num_vel_pts_blade) + 1);
                        view.vel_pos = ::kynema_sgf::utils::slice(
                            vel_pos_now, 0, tdata.num_vel_pts_blade);
                    }
                    // Second input value of epsilon_chord, in the solver frame
                    constexpr int ic = epsilon_component(
                        ::ext_turb::blade_frame_of<SolverTurbine>(), 1);
                    fllc_init(
                        tdata.fllc[i], view, tdata.eps_chord[ic],
                        data.info().label + " blade " + std::to_string(i),
                        data.info().is_root_proc);
                    started = true;
                }
            }

            if (started && restart_call && data.info().is_root_proc) {
                if (!external_read_checkpoint<DataType>(
                        data, io.restart_file())) {
                    amrex::OutStream()
                        << "\nWARNING: Actuator " << data.info().label
                        << ": the state of the filtered lifting line "
                           "correction "
                        << "in the checkpoint " << io.restart_file()
                        << " does not fit this turbine: it is for other "
                           "numbers "
                        << "of blades or points, for other positions of the "
                        << "points or for other inputs of the correction. The "
                        << "correction starts from zero." << std::endl;
                }
            }
            if (waits && state_in_checkpoint && data.info().is_root_proc) {
                amrex::OutStream()
                    << "\nWARNING: Actuator " << data.info().label
                    << ": the external solver starts this turbine anew, so the "
                    << "state of the filtered lifting line correction in the "
                    << "checkpoint " << io.restart_file()
                    << " is not used. The correction starts from zero in the "
                    << "first time step." << std::endl;
            } else if (waits && data.info().is_root_proc) {
                amrex::OutStream()
                    << "\nWARNING: Actuator " << data.info().label
                    << ": the checkpoint " << io.restart_file()
                    << " has no state of the filtered lifting line correction "
                    << "for this time. If the correction had not started when "
                    << "the checkpoint was written, nothing is lost. If it had "
                    << "(a checkpoint written without the state), the restart "
                    << "loses it. The correction starts from zero in the first "
                    << "time step." << std::endl;
            }
        }
        tdata.at_first_force_call = false;
    }
};

} // namespace kynema_sgf::actuator::ops

#endif /* TURBINE_EXTERNAL_OPS_H */
