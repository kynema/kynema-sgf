#ifndef TURBINE_EXTERNAL_OPS_H
#define TURBINE_EXTERNAL_OPS_H

#include "amr-wind/wind_energy/actuator/turbine/turbine_utils.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/wind_energy/actuator/actuator_utils.H"
#include "amr-wind/wind_energy/actuator/FLLCOp.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/wind_energy/actuator/turbine/external/turbine_external_utils.H"

using namespace amr_wind::actuator::external;

namespace amr_wind::actuator::ops {

template <typename datatype>
void read_ops(datatype& data, const utils::ActParser& pp)
{
    // Tell reading utility the name of the solver
    data.meta().solver_name = data.meta().ext_ptr->identifier();
    // Data common to any turbine actuator simulation
    utils::read_inputs(data.meta(), data.info(), pp);

    auto& tdata = data.meta();

    // Get density value for normalization
    pp.get("density", tdata.density);

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

    const auto& time = data.sim().time();
    tf.chkpt_interval = time.chkpt_interval();
}

template <typename datatype>
void external_determine_influenced_procs(datatype& data)
{
    auto& info = data.info();
    info.procs =
        utils::determine_influenced_procs(data.sim().mesh(), info.bound_box);

    AMREX_ALWAYS_ASSERT(info.root_proc > -1);
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
    namespace utils = ::amr_wind::actuator::utils;
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

    const auto& tdata = data.meta();
    const auto& bp = data.info().base_pos;
    const auto& pxvel = tdata.ext_data.position_at_vel(0);
    const auto& pyvel = tdata.ext_data.position_at_vel(1);
    const auto& pzvel = tdata.ext_data.position_at_vel(2);
    auto& vel_pos = data.grid().vel_pos;
    for (int i = 0; i < vel_pos.size(); ++i) {
        vel_pos[i].x() = static_cast<amrex::Real>(pxvel[i]) + bp.x();
        vel_pos[i].y() = static_cast<amrex::Real>(pyvel[i]) + bp.y();
        vel_pos[i].z() = static_cast<amrex::Real>(pzvel[i]) + bp.z();
    }
}

template <typename datatype>
void update_vel_op(datatype& data)
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
        // Compute the relative velocity needed for the FLLC
        const auto& xdot = tdata.ext_data.solid_velocity(0);
        const auto& ydot = tdata.ext_data.solid_velocity(1);
        const auto& zdot = tdata.ext_data.solid_velocity(2);
        for (int i = 0; i < tdata.vel_rel.size(); ++i) {
            tdata.vel_rel[i][0] = uvel[i] - static_cast<amrex::Real>(xdot[i]);
            tdata.vel_rel[i][1] = vvel[i] - static_cast<amrex::Real>(ydot[i]);
            tdata.vel_rel[i][2] = wvel[i] - static_cast<amrex::Real>(zdot[i]);
        }
        // Loop through each blade and apply the FLLC
        for (int i = 0; i < tdata.num_blades; ++i) {
            FLLCOp()(tdata.blades[i], tdata.fllc[i]);
        }
    }

    for (int i = 0; i < tdata.ext_data.length_fluid_velocity(0); ++i) {
        uvel[i] = static_cast<float>(vel[i].x());
        vvel[i] = static_cast<float>(vel[i].y());
        wvel[i] = static_cast<float>(vel[i].z());
    }
}

template <typename datatype>
void compute_force_op(datatype& data)
{
    // Advance external solver by specified number of sub-steps
    ext_step<datatype>(data);
    // Broadcast data to all the processes that contain patches influenced
    // by this turbine
    scatter_data<datatype>(data);

    const auto& time = data.sim().time();

    auto& tdata = data.meta();
    if (!tdata.fllc.empty()) {
        for (int i = 0; i < tdata.num_blades; ++i) {
            if (!(tdata.fllc[i].initialized) &&
                (time.current_time() > tdata.fllc[i].fllc_start_time)) {
                fllc_init(tdata.fllc[i], tdata.blades[i], tdata.eps_chord[0]);
            }
        }
    }
}

} // namespace amr_wind::actuator::ops

#endif /* TURBINE_EXTERNAL_OPS_H */
