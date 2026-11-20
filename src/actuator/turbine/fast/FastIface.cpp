#include "src/actuator/turbine/external/ExtTurbIface.H"
#include "src/actuator/turbine/fast/fast_types.H"
#include "src/actuator/turbine/fast/fast_wrapper.H"
#include "src/CFDSim.H"
#include "src/core/SimTime.H"
#include "src/utilities/io_utils.H"
#include "src/utilities/ncutils/nc_interface.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_FileSystem.H"

#include <algorithm>
#include <cmath>
#include "AMReX_REAL.H"

using namespace amrex::literals;

using namespace sgf_fast;

namespace ext_turb {
namespace {

template <typename FType, class... Args>
void fast_func(const FType&& func, Args... args)
{
    int ierr = ErrID_None;
    amrex::Array<char, fast_strlen()> err_msg;
    func(std::forward<Args>(args)..., &ierr, err_msg.begin());
    if (ierr >= ErrID_Fatal) {
        std::string prefix = "FastIface: Error calling OpenFAST function: \n";
        amrex::Abort(prefix + std::string(err_msg.begin()));
    }
}

void copy_filename(const std::string& inp, char* out)
{
    const int str_len = static_cast<int>(inp.size());
    if (str_len >= fast_strlen()) {
        amrex::Abort(
            "FastIface: Filename greater than maximum allowable length: " +
            inp);
    }
    const int len = amrex::min(fast_strlen() - 1, str_len);
    std::copy(inp.data(), inp.data() + len, out);
    out[len] = '\0';
}

} // namespace

template <>
ExtTurbIface<FastTurbine, FastSolverData>::~ExtTurbIface()
{
    int ierr = ErrID_None;
    amrex::Array<char, fast_strlen()> err_msg;
    FAST_DeallocateTurbines(&ierr, err_msg.begin());

    if (ierr != ErrID_None) {
        amrex::OutStream() << "\nFastIface: Error deallocating turbine data\n"
                           << err_msg.begin() << std::endl;
    }
}

template <>
void ExtTurbIface<FastTurbine, FastSolverData>::parse_inputs(
    const kynema_sgf::CFDSim& sim, const std::string& inp_name)
{
    amrex::ParmParse pp(inp_name);

    const auto& time = sim.time();
    // Check that the user has not enabled adaptive timestepping
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !time.adaptive_timestep(),
        "Adaptive time-stepping not supported when OpenFAST is enabled");

    m_dt_cfd = time.delta_t();

    // Set OpenFAST end time to be at least as long as the CFD time. User
    // can choose a longer duration in input file.
    const amrex::Real stop1 = time.stop_time() > 0.0_rt
                                  ? time.stop_time()
                                  : std::numeric_limits<amrex::Real>::max();
    const amrex::Real stop2 = time.stop_time_index() > 0
                                  ? time.stop_time_index() * time.delta_t()
                                  : std::numeric_limits<amrex::Real>::max();
    const amrex::Real cfd_stop = amrex::min<amrex::Real>(stop1, stop2);
    m_stop_time = cfd_stop;

    pp.query("start_time", m_start_time);
    pp.query("stop_time", m_stop_time);

    // Ensure that the user specified m_stop_time is not shorter than CFD sim
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_stop_time > (cfd_stop - 1.0e-6_rt),
        "OpenFAST stop_time = " + std::to_string(m_stop_time) +
            ", CFD stop_time = " + std::to_string(cfd_stop));

    if (m_start_time > 0.0_rt) {
        m_sim_mode = ::ext_turb::SimMode::replay;

        std::string sim_mode{"replay"};
        pp.query("openfast_sim_mode", sim_mode);
        if (sim_mode == "replay") {
            m_sim_mode = ::ext_turb::SimMode::replay;
        } else if (sim_mode == "restart") {
            m_sim_mode = ::ext_turb::SimMode::restart;
        } else {
            amrex::Abort(
                "Invalid simulation mode when start time > 0 provided: " +
                sim_mode);
        }
    }
}

template <>
void ExtTurbIface<FastTurbine, FastSolverData>::allocate_ext_turbines()
{
    BL_PROFILE("kynema-sgf::FastIface::allocate_turbines");
    // OpenFAST keeps one count of time steps for all the turbines of a
    // process and moves it on when the last of them is stepped
    // (FAST_CFD_Step), while a turbine here takes all its steps of a CFD
    // time step, and writes its checkpoint, before the next turbine is
    // stepped; a replay takes all the steps of one turbine in a row. With
    // several turbines in a process, a turbine that takes more than one
    // OpenFAST step per CFD step repeats a time step or gets the time of
    // another turbine, and the checkpoint of every turbine but the last is
    // named and stamped one step low. So every OpenFAST turbine needs a
    // process of its own. Every turbine of this process is registered by now
    // (register_turbine refuses one after the allocation), so the run stops
    // here, before OpenFAST reads its inputs.
    if (m_turbine_data.size() > 1) {
        amrex::Abort(
            "FastIface: process " +
            std::to_string(amrex::ParallelDescriptor::MyProc()) + " holds " +
            std::to_string(m_turbine_data.size()) +
            " OpenFAST turbines. OpenFAST counts the time steps of the "
            "turbines of a process together, so they cannot be stepped, "
            "checkpointed or replayed one after the other: use at least as "
            "many MPI ranks as there are actuators, so that every turbine has "
            "a rank of its own.");
    }
    int nturbines = static_cast<int>(m_turbine_data.size());
    fast_func(FAST_AllocateTurbines, &nturbines);
    m_is_initialized = true;
}

template <>
void ExtTurbIface<FastTurbine, FastSolverData>::init_solution(
    const int local_id)
{
    BL_PROFILE("kynema-sgf::FastIface::init_solution");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        local_id < static_cast<int>(m_turbine_data.size()),
        "local_id = " + std::to_string(local_id) +
            ", m_turbine_data.size() = " +
            std::to_string(m_turbine_data.size()));
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_is_initialized,
        "m_is_initialized = " + std::to_string(m_is_initialized));

    auto& fi = *m_turbine_data[local_id];
    fast_func(FAST_Solution0, &fi.tid_local);
    fi.is_solution0 = false;
}

template <>
void ExtTurbIface<FastTurbine, FastSolverData>::get_hub_stats(
    const int local_id)
{
    BL_PROFILE("kynema-sgf::FastIface::get_hub_stats");

    auto& fi = *m_turbine_data[local_id];

    fast_func(
        FAST_HubPosition, &fi.tid_local, fi.hub_abs_pos.begin(),
        fi.hub_rot_vel.begin(), fi.hub_orient.begin());
}

#ifdef KYNEMA_SGF_USE_OPENFAST

template <>
void ExtTurbIface<FastTurbine, FastSolverData>::prepare_netcdf_file(
    FastTurbine& fi)
{
#ifdef KYNEMA_SGF_USE_NETCDF
    BL_PROFILE("kynema-sgf::FastIface::prepare_netcdf_file");
    if (!amrex::UtilCreateDirectory(m_solver_data.output_dir, 0755)) {
        amrex::CreateDirectoryFailed(m_solver_data.output_dir);
    }

    const std::string fname =
        m_solver_data.output_dir + "/" + fi.tlabel + ".nc";

    // Don't overwrite existing
    if (amrex::FileSystem::Exists(fname)) {
        return;
    }

    auto ncf = ncutils::NCFile::create(fname, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    const std::string np_name = "num_vel_points";
    ncf.enter_def_mode();
    ncf.put_attr("title", "Kynema-SGF OpenFAST velocity data");
    ncf.put_attr(
        "Kynema-SGF_version", kynema_sgf::ioutils::kynema_sgf_version());
    ncf.put_attr("created_on", kynema_sgf::ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim(np_name, fi.from_cfd.u_Len);
    ncf.def_dim("ndim", AMREX_SPACEDIM);
    ncf.def_var("time", NC_FLOAT, {nt_name});
    ncf.def_var("xco", NC_FLOAT, {np_name});
    ncf.def_var("yco", NC_FLOAT, {np_name});
    ncf.def_var("zco", NC_FLOAT, {np_name});
    ncf.def_var("uvel", NC_FLOAT, {nt_name, np_name});
    ncf.def_var("vvel", NC_FLOAT, {nt_name, np_name});
    ncf.def_var("wvel", NC_FLOAT, {nt_name, np_name});
    ncf.exit_def_mode();

    {
        const auto npts = static_cast<size_t>(fi.from_cfd.u_Len);
        auto xco = ncf.var("xco");
        xco.put(fi.position_at_vel(0), {0}, {npts});
        auto yco = ncf.var("yco");
        yco.put(fi.position_at_vel(1), {0}, {npts});
        auto zco = ncf.var("zco");
        zco.put(fi.position_at_vel(2), {0}, {npts});
    }
#else
    amrex::ignore_unused(fi);
    amrex::OutStream()
        << "WARNING: FastIface: NetCDF support was not enabled during compile "
           "time. FastIface cannot support legacy replay capability."
        << std::endl;
#endif
}

template <>
void ExtTurbIface<FastTurbine, FastSolverData>::write_velocity_data(
    const FastTurbine& fi)
{
#ifdef KYNEMA_SGF_USE_NETCDF
    BL_PROFILE("kynema-sgf::FastIface::write_velocity_data");
    const std::string fname =
        m_solver_data.output_dir + "/" + fi.tlabel + ".nc";
    auto ncf = ncutils::NCFile::open(fname, NC_WRITE);
    // Record n holds the velocities of CFD time step n of the turbine,
    // counted from the start of the external solver, which is how a replay
    // reads the file. A run restarted from a checkpoint older than the last
    // record writes over the records of the steps it takes again.
    const auto nt = static_cast<size_t>(fi.time_index / fi.num_substeps);
    const auto npts = static_cast<size_t>(fi.from_cfd.u_Len);

    const amrex::Real time = fi.time_index * fi.dt_ext;
    ncf.var("time").put(&time, {nt}, {1});
    const auto& uu = fi.from_cfd;
    ncf.var("uvel").put(uu.u, {nt, 0}, {1, npts});
    ncf.var("vvel").put(uu.v, {nt, 0}, {1, npts});
    ncf.var("wvel").put(uu.w, {nt, 0}, {1, npts});
#else
    amrex::ignore_unused(fi);
#endif
}

template <>
void ExtTurbIface<FastTurbine, FastSolverData>::read_velocity_data(
    FastTurbine& fi, const ncutils::NCFile& ncf, const size_t tid)
{
#ifdef KYNEMA_SGF_USE_NETCDF
    const auto nt = static_cast<size_t>(tid);
    const auto npts = static_cast<size_t>(fi.from_cfd.u_Len);

    // Record nt holds the velocities of CFD time step nt of the turbine. A
    // record of another time, or one never written, would be replayed without
    // notice.
    amrex::Real time = 0.0_rt;
    ncf.var("time").get(&time, {nt}, {1});
    const amrex::Real dt_record =
        static_cast<amrex::Real>(fi.num_substeps) * fi.dt_ext;
    const amrex::Real time_expected = static_cast<amrex::Real>(nt) * dt_record;
    const bool never_written =
        !(time < 0.5_rt * static_cast<amrex::Real>(NC_FILL_FLOAT));
    if (never_written || !(std::abs(time - time_expected) <=
                           0.5_rt * dt_record + 1.0e-6_rt * time_expected)) {
        amrex::Abort(
            "FastIface: record " + std::to_string(nt) +
            " of the velocity file of " + fi.tlabel +
            (never_written ? " was never written"
                           : " has time " + std::to_string(time)) +
            ", time expected " + std::to_string(time_expected) +
            ". The file is not the one written by the run being replayed");
    }

    auto& uu = fi.from_cfd;
    ncf.var("uvel").get(uu.u, {nt, 0}, {1, npts});
    ncf.var("vvel").get(uu.v, {nt, 0}, {1, npts});
    ncf.var("wvel").get(uu.w, {nt, 0}, {1, npts});
#else
    amrex::ignore_unused(fi);
    amrex::Abort(
        "FastIface::read_velocity_data: Kynema-SGF was not compiled with "
        "NetCDF "
        "support");
#endif
}

#else

template <>
void ExtTurbIface<FastTurbine, FastSolverData>::prepare_netcdf_file(
    FastTurbine& /*unused*/)
{}
template <>
void ExtTurbIface<FastTurbine, FastSolverData>::write_velocity_data(
    const FastTurbine& /*unused*/)
{}
template <>
void ExtTurbIface<FastTurbine, FastSolverData>::read_velocity_data(
    FastTurbine& /*unused*/,
    const ncutils::NCFile& /*unused*/,
    const size_t /*unused*/)
{}

#endif

template <>
void ExtTurbIface<FastTurbine, FastSolverData>::calc_aero_loads(
    FastTurbine& /*fi*/)
{
    // Nothing for OpenFAST as aerodynamic loads cannot be calculated
    // separately.
}

template <>
void ExtTurbIface<FastTurbine, FastSolverData>::do_turbine_step(FastTurbine& fi)
{
    // One OpenFAST turbine per process (allocate_ext_turbines)
    AMREX_ALWAYS_ASSERT(m_turbine_data.size() == 1);
    fast_func(FAST_Step, &fi.tid_local);
}

template <>
void ExtTurbIface<FastTurbine, FastSolverData>::write_turbine_checkpoint(
    int& tid)
{
    amrex::Array<char, fast_strlen()> rst_file;
    copy_filename(" ", rst_file.begin());
    fast_func(FAST_CreateCheckpoint, &tid, rst_file.begin());
}

template <>
void ExtTurbIface<FastTurbine, FastSolverData>::ext_init_turbine(
    FastTurbine& fi)
{
    BL_PROFILE("kynema-sgf::FastIface::init_turbine");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        amrex::FileSystem::Exists(fi.input_file),
        "FastIface: Cannot find OpenFAST input file: " + fi.input_file);

    int abort_lev;
    amrex::Array<char, fast_strlen()> inp_file;
    copy_filename(fi.input_file, inp_file.begin());

#if OPENFAST_VERSION_MAJOR >= 4
    char out_file[fast_strlen()];
    fast_func(
        FAST_ExtInfw_Init, &fi.tid_local, &fi.stop_time, inp_file.begin(),
        &fi.tid_global, out_file, &fi.num_pts_blade, &fi.num_pts_tower,
        fi.base_pos.begin(), &abort_lev, &fi.dt_cfd, &fi.dt_ext,
        &m_solver_data.m_inflow_type, &fi.num_blades, &fi.num_blade_elem,
        &fi.num_tower_elem, &fi.chord_cluster_type, &fi.to_cfd, &fi.from_cfd);
#else
    fast_func(
        FAST_OpFM_Init, &fi.tid_local, &fi.stop_time, inp_file.begin(),
        &fi.tid_global, &m_solver_data.m_num_sc_inputs_glob,
        &m_solver_data.m_num_sc_inputs, &m_solver_data.m_num_sc_outputs,
        &m_solver_data.m_init_sc_inputs_glob,
        &m_solver_data.m_init_sc_inputs_turbine, &fi.num_pts_blade,
        &fi.num_pts_tower, fi.base_pos.begin(), &abort_lev, &fi.dt_ext,
        &fi.num_blades, &fi.num_blade_elem, &fi.chord_cluster_type, &fi.to_cfd,
        &fi.from_cfd, &fi.to_sc, &fi.from_sc);
#endif

    {
#ifdef KYNEMA_SGF_USE_OPENFAST
        // Check if OpenFAST has tower and reset tower nodes appropriately
        const int npts = fi.to_cfd.fx_Len;
        const int nrotor_pts = fi.num_blades * fi.num_pts_blade + 1;
        if (nrotor_pts == npts) {
            amrex::OutStream()
                << "OpenFAST model does not include tower for turbine: "
                << fi.tlabel << " Turning off tower actuator points"
                << std::endl;
            fi.num_pts_tower = 0;
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            npts == (nrotor_pts + fi.num_pts_tower),
            "OpenFAST points: npts = " + std::to_string(npts) +
                ", nrotor_pts = " + std::to_string(nrotor_pts) +
                ", num_pts_tower = " + std::to_string(fi.num_pts_tower));
#endif
    }

    // OpenFAST steps per CFD step
    fi.num_substeps = num_whole_steps(
        fi.dt_cfd, fi.dt_ext, "FastIface " + fi.tlabel + ": time.fixed_dt",
        "the OpenFAST DT");
    AMREX_ALWAYS_ASSERT(fi.num_substeps > 0);
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
template <>
void ExtTurbIface<FastTurbine, FastSolverData>::ext_replay_turbine(
    FastTurbine& fi)
{
#ifdef KYNEMA_SGF_USE_NETCDF
    BL_PROFILE("kynema-sgf::FastIface::replay_turbine");

    // CFD steps up to the start time, one velocity record each. Tolerance
    // of 1e-3 step for a start time copied from an accumulated CFD time.
    const int num_cfd_steps = num_whole_steps(
        fi.start_time, fi.dt_cfd,
        "FastIface " + fi.tlabel + ": Actuator." + fi.tlabel +
            ".openfast_start_time",
        "time.fixed_dt", 1.0e-3_rt);
    // Number of OpenFAST time steps this turbine is advanced by
    const int num_steps = num_cfd_steps * fi.num_substeps;

    // Ensure that the NetCDF file exists and contains the required number of
    // timesteps for restart.
    const std::string fname =
        m_solver_data.output_dir + "/" + fi.tlabel + ".nc";
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        amrex::FileSystem::Exists(fname),
        "FastIface: Cannot find OpenFAST velocity data file: " + fname);

    auto ncf = ncutils::NCFile::open(fname, NC_NOWRITE);
    {
        const int nt = static_cast<int>(ncf.dim("num_time_steps").len());
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            nt >= num_cfd_steps,
            "OpenFAST velocity timesteps = " + std::to_string(nt) +
                ", required CFD timesteps = " + std::to_string(num_cfd_steps));
    }

    // The replay takes all the time steps of this turbine in a row: one
    // OpenFAST turbine per process (allocate_ext_turbines)
    AMREX_ALWAYS_ASSERT(m_turbine_data.size() == 1);

    // Replay OpenFAST simulation for the desired number of timesteps to mimic
    // restart
    fi.time_index = 0;
    read_velocity_data(fi, ncf, 0);
    fast_func(FAST_Solution0, &fi.tid_local);
    fi.is_solution0 = false;

    for (int ic = 0; ic < num_cfd_steps; ++ic) {
        read_velocity_data(fi, ncf, ic);
        for (int i = 0; i < fi.num_substeps; ++i, ++fi.time_index) {
            fast_func(FAST_Step, &fi.tid_local);
        }
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        fi.time_index == num_steps,
        "fi.time_index = " + std::to_string(fi.time_index) +
            ", num_steps = " + std::to_string(num_steps));
#else
    amrex::ignore_unused(fi);
    amrex::Abort(
        "FastIface::replay_turbine: Kynema-SGF was not compiled with NetCDF "
        "support");
#endif
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
template <>
void ExtTurbIface<FastTurbine, FastSolverData>::ext_restart_turbine(
    FastTurbine& fi)
{
    BL_PROFILE("kynema-sgf::FastIface::restart_turbine");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        amrex::FileSystem::Exists(fi.checkpoint_file + ".chkp"),
        "FastIface: Cannot find OpenFAST checkpoint file: " +
            fi.checkpoint_file);

    int abort_lev;
    amrex::Array<char, fast_strlen()> chkpt_file;
    copy_filename(fi.checkpoint_file, chkpt_file.begin());

#if OPENFAST_VERSION_MAJOR >= 4
    fast_func(
        FAST_ExtInfw_Restart, &fi.tid_local, chkpt_file.begin(), &abort_lev,
        &fi.dt_ext, &fi.num_blades, &fi.num_blade_elem, &fi.num_tower_elem,
        &fi.time_index, &fi.to_cfd, &fi.from_cfd);
#else
    fast_func(
        FAST_OpFM_Restart, &fi.tid_local, chkpt_file.begin(), &abort_lev,
        &fi.dt_ext, &fi.num_blades, &fi.num_blade_elem, &fi.time_index,
        &fi.to_cfd, &fi.from_cfd, &fi.to_sc, &fi.from_sc);
#endif

    {
#ifdef KYNEMA_SGF_USE_OPENFAST
        // Check if OpenFAST has tower and reset tower nodes appropriately
        const int npts = fi.to_cfd.fx_Len;
        const int nrotor_pts = fi.num_blades * fi.num_pts_blade + 1;
        if (nrotor_pts == npts) {
            amrex::OutStream()
                << "OpenFAST model does not include tower for turbine: "
                << fi.tlabel << " Turning off tower actuator points"
                << std::endl;
            fi.num_pts_tower = 0;
        }
        // OpenFAST takes the number of actuator points from its checkpoint,
        // not from the inputs of this run
        if (npts != (nrotor_pts + fi.num_pts_tower)) {
            amrex::Abort(
                "FastIface: the OpenFAST checkpoint of " + fi.tlabel + " has " +
                std::to_string(npts) +
                " actuator points (blades, hub and tower); this run asks for " +
                std::to_string(nrotor_pts + fi.num_pts_tower) +
                " (num_points_blade = " + std::to_string(fi.num_pts_blade) +
                ", num_points_tower = " + std::to_string(fi.num_pts_tower) +
                "). A restart keeps the number of actuator points of the run "
                "that wrote the checkpoint.");
        }
#endif
    }

    // OpenFAST steps per CFD step
    fi.num_substeps = num_whole_steps(
        fi.dt_cfd, fi.dt_ext, "FastIface " + fi.tlabel + ": time.fixed_dt",
        "the OpenFAST DT");
    AMREX_ALWAYS_ASSERT(fi.num_substeps > 0);
}

template class ExtTurbIface<FastTurbine, FastSolverData>;

} // namespace ext_turb
