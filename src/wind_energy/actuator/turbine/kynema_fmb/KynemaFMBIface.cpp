#include <algorithm>
#include <cmath>
#include "src/wind_energy/actuator/turbine/external/ExtTurbIface.H"
#include "src/wind_energy/actuator/turbine/kynema_fmb/kynema_fmb_types.H"
#include "src/CFDSim.H"
#include "src/core/SimTime.H"
#include "src/utilities/io_utils.H"
#include "src/utilities/ncutils/nc_interface.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_FileSystem.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace ext_turb {

template <>
ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::~ExtTurbIface()
{
    //! Do deallocation if necessary
}

// !! This doesn't get used by the actual code, just by unit tests !! //
template <>
void ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::parse_inputs(
    const kynema_sgf::CFDSim& sim, const std::string& inp_name)
{
    amrex::ParmParse pp(inp_name);

    const auto& time = sim.time();
    //! Check that the user has not enabled adaptive timestepping
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !time.adaptive_timestep(),
        "Adaptive time-stepping not supported when Kynema is enabled");

    m_dt_cfd = time.delta_t();
}

template <>
void ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::
    allocate_ext_turbines()
{
    BL_PROFILE("kynema-sgf::KynemaFMBIface::allocate_turbines");

    m_is_initialized = true;

    // Get solver parameters
    {
        amrex::ParmParse pp("KynemaFMB");
        pp.query("damping_factor", m_solver_data.damping_factor);
        pp.query("max_nonlinear_iterations", m_solver_data.nl_iter_max);
        pp.query("abs_err_tol", m_solver_data.abs_err_tol);
        pp.query("rel_err_tol", m_solver_data.rel_err_tol);
    }

    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_solver_data.gravity);
    }
}

template <>
void ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::init_solution(
    const int local_id)
{
    BL_PROFILE("kynema-sgf::KynemaFMBIface::init_solution");
    AMREX_ALWAYS_ASSERT(local_id < static_cast<int>(m_turbine_data.size()));
    AMREX_ALWAYS_ASSERT(m_is_initialized);

    auto& fi = *m_turbine_data[local_id];

    ::sgf_fmb::update_turbine_model(fi, false);

    fi.is_solution0 = false;
}

template <>
void ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::get_hub_stats(
    const int local_id)
{
    BL_PROFILE("kynema-sgf::KynemaFMBIface::get_hub_stats");

    auto& fi = *m_turbine_data[local_id];

    ::sgf_fmb::get_hub_state(fi);
    // Orientation of hub is already in turbine data as first point
    for (int comp = 0; comp < AMREX_SPACEDIM * AMREX_SPACEDIM; ++comp) {
        fi.hub_orient[comp] = fi.orientation()[comp];
    }
}

#ifdef KYNEMA_SGF_USE_KYNEMA_FMB

template <>
void ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::prepare_netcdf_file(
    KynemaFMBTurbine& fi)
{
#ifdef KYNEMA_SGF_USE_NETCDF
    BL_PROFILE("kynema-sgf::KynemaFMBIface::prepare_netcdf_file");
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
    ncf.put_attr("title", "Kynema-SGF Kynema velocity data");
    ncf.put_attr(
        "Kynema-SGF_version", kynema_sgf::ioutils::kynema_sgf_version());
    ncf.put_attr("created_on", kynema_sgf::ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim(np_name, fi.length_fluid_velocity(0));
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
        const auto npts = static_cast<size_t>(fi.length_fluid_velocity(0));
        auto xco = ncf.var("xco");
        xco.put(fi.position_at_vel(0), {0}, {npts});
        auto yco = ncf.var("yco");
        yco.put(fi.position_at_vel(1), {0}, {npts});
        auto zco = ncf.var("zco");
        zco.put(fi.position_at_vel(2), {0}, {npts});
    }
#else
    amrex::ignore_unused(fi);
    amrex::OutStream() << "WARNING: KynemaFMBIface: NetCDF support was not "
                          "enabled during compile "
                          "time. KynemaFMBIface cannot support restart."
                       << std::endl;
#endif
}

template <>
void ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::write_velocity_data(
    const KynemaFMBTurbine& fi)
{
#ifdef KYNEMA_SGF_USE_NETCDF
    BL_PROFILE("kynema-sgf::KynemaFMBIface::write_velocity_data");
    const std::string fname =
        m_solver_data.output_dir + "/" + fi.tlabel + ".nc";
    auto ncf = ncutils::NCFile::open(fname, NC_WRITE);
    const std::string nt_name = "num_time_steps";
    const size_t nt = ncf.dim(nt_name).len();
    const auto npts = static_cast<size_t>(fi.length_fluid_velocity(0));

    const amrex::Real time = fi.time_index * fi.dt_ext;
    ncf.var("time").put(&time, {nt}, {1});
    ncf.var("uvel").put(fi.fluid_velocity(0), {nt, 0}, {1, npts});
    ncf.var("vvel").put(fi.fluid_velocity(1), {nt, 0}, {1, npts});
    ncf.var("wvel").put(fi.fluid_velocity(2), {nt, 0}, {1, npts});
#else
    amrex::ignore_unused(fi);
#endif
}

template <>
void ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::read_velocity_data(
    KynemaFMBTurbine& fi, const ncutils::NCFile& ncf, const size_t tid)
{
#ifdef KYNEMA_SGF_USE_NETCDF
    const auto nt = static_cast<size_t>(tid);
    const auto npts = static_cast<size_t>(fi.length_fluid_velocity(0));

    ncf.var("uvel").get(fi.fluid_velocity(0), {nt, 0}, {1, npts});
    ncf.var("vvel").get(fi.fluid_velocity(1), {nt, 0}, {1, npts});
    ncf.var("wvel").get(fi.fluid_velocity(2), {nt, 0}, {1, npts});
#else
    amrex::ignore_unused(fi);
    amrex::Abort(
        "KynemaFMBIface::read_velocity_data: Kynema-SGF was not compiled with "
        "NetCDF "
        "support");
#endif
}

#else

template <>
void ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::prepare_netcdf_file(
    KynemaFMBTurbine& /*unused*/)
{}
template <>
void ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::write_velocity_data(
    const KynemaFMBTurbine& /*unused*/)
{}
template <>
void ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::read_velocity_data(
    KynemaFMBTurbine& /*unused*/,
    const ncutils::NCFile& /*unused*/,
    const size_t /*unused*/)
{}

#endif

template <>
void ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::calc_aero_loads(
    KynemaFMBTurbine& fi)
{
    ::sgf_fmb::calculate_aerodynamic_loads(fi);
}

template <>
void ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::do_turbine_step(
    KynemaFMBTurbine& fi)
{
    ::sgf_fmb::update_turbine_model(fi, true);
}

template <>
void ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::
    write_turbine_checkpoint(int& tid)
{
    BL_PROFILE("kynema-sgf::KynemaFMBIface::write_turbine_checkpoint");
    auto& fi = *m_turbine_data[tid];

    // Write the checkpoint file
    ::sgf_fmb::write_checkpoint_file(
        fi, sgf_fmb::checkpoint_filename(
                fi.time_index / fi.num_substeps, fi.tlabel));
}

template <>
void ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::ext_init_turbine(
    KynemaFMBTurbine& fi)
{
    BL_PROFILE("kynema-sgf::KynemaFMBIface::init_turbine");

    ::sgf_fmb::init_turbine_interface(fi, m_solver_data);

    // Determine the number of substeps for Kynema per CFD timestep
    fi.num_substeps = static_cast<int>(std::floor(fi.dt_cfd / fi.dt_ext));

    AMREX_ALWAYS_ASSERT(fi.num_substeps > 0);
    // Check that the time step sizes are consistent and Kynema advances at an
    // integral multiple of CFD timestep
    amrex::Real dt_err =
        fi.dt_cfd / (static_cast<amrex::Real>(fi.num_substeps) * fi.dt_ext) -
        1.0_rt;
    if (dt_err > 1.0e-12_rt) {
        amrex::Abort(
            "KynemaFMBIFace: Kynema-FMB timestep is not an integral "
            "multiple of CFD timestep");
    }

    // Offset time index if restarting from a checkpoint file
    if (!fi.checkpoint_file.empty()) {
        fi.time_index *= fi.num_substeps;
        ::sgf_fmb::set_time_step_index(fi);
    }

    fi.allocate_buffers();
    ::sgf_fmb::update_turbine_model(fi, false);
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
template <>
void ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::ext_replay_turbine(
    KynemaFMBTurbine& /*fi*/)
{
    // Code should never reach here
    amrex::Abort(
        "KynemaFMBIface: replaying Kynema-FMB turbine is not implemented");
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
template <>
void ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>::ext_restart_turbine(
    KynemaFMBTurbine& fi)
{
    BL_PROFILE("kynema-sgf::KynemaFMBIface::restart_turbine");

    ext_init_turbine(fi);
}

template class ExtTurbIface<KynemaFMBTurbine, KynemaFMBSolverData>;

} // namespace ext_turb
