#include "amr-wind/wind_energy/actuator/turbine/turbine_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/wind_energy/actuator/FLLC.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator::utils {

void read_inputs(
    TurbineBaseData& tdata, TurbineInfo& tinfo, const utils::ActParser& pp)
{
    pp.query("num_blades", tdata.num_blades);
    pp.get("num_points_blade", tdata.num_pts_blade);
    tdata.num_vel_pts_blade = tdata.num_pts_blade;
    pp.get("num_points_tower", tdata.num_pts_tower);
    pp.query("nacelle_area", tdata.nacelle_area);
    pp.query("nacelle_drag_coeff", tdata.nacelle_cd);

    if (!pp.contains("epsilon") && !pp.contains("epsilon_chord")) {
        amrex::Abort(
            "Actuator turbine simulations require specification of one or both "
            "of 'epsilon' or 'epsilon_chord'");
    }

    pp.query_either("epsilon", tdata.eps_inp);
    pp.query_either("epsilon_chord", tdata.eps_chord);
    pp.query_either("epsilon_min", tdata.eps_min);
    pp.query_either("epsilon_tower", tdata.eps_tower);

    pp.get("base_position", tinfo.base_pos);
    if (tdata.solver_name != "Kynema") {
        pp.get("rotor_diameter", tinfo.rotor_diameter);
        pp.get("hub_height", tinfo.hub_height);
    }
    bool use_fllc = false;
    pp.query("fllc", use_fllc);
    if (use_fllc) {
        std::string type;
        pp.query("type", type);
        if (type == "TurbineFastDisk") {
            amrex::Abort(
                "FLLC is not compatible with TurbineFastDisk turbine type");
        }
        for (int i = 0; i < tdata.num_blades; ++i) {
            tdata.fllc.emplace_back();
            fllc_parse(pp, tdata.fllc.back());
        }
    }

    // clang-format off
    const auto& bp = tinfo.base_pos;
    const auto& rad = 0.5_rt * tinfo.rotor_diameter;
    const auto& hh = tinfo.hub_height;
    tinfo.bound_box = amrex::RealBox(
        bp.x() - (1.25_rt * rad), bp.y() - (1.25_rt * rad), bp.z() - (1.25_rt * rad),
        bp.x() + (1.25_rt * rad), bp.y() + (1.25_rt * rad), bp.z() + (1.25_rt * rad) + hh
    );
    // clang-format on
}

void prepare_netcdf_file(
    const std::string& ncfile,
    const TurbineBaseData& meta,
    const TurbineInfo& info,
    const ActGrid& grid)
{
#ifdef AMR_WIND_USE_NETCDF
    // Only root process handles I/O
    if (!info.is_root_proc) {
        return;
    }

    auto ncf = ncutils::NCFile::create(ncfile, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    const std::string np_name = "num_actuator_points";
    const std::string nvel_name = "num_vel_points";
    const std::vector<std::string> two_dim{nt_name, np_name};

    ncf.enter_def_mode();
    ncf.put_attr("title", "AMR-Wind turbine actuator output");
    ncf.put_attr("version", ioutils::amr_wind_version());
    ncf.put_attr("created_on", ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim("ndim", AMREX_SPACEDIM);
    ncf.def_dim(
        "mat_dim", static_cast<size_t>(AMREX_SPACEDIM * AMREX_SPACEDIM));

    auto grp = ncf.def_group(info.label);
    grp.put_attr("num_blades", std::vector<int>{meta.num_blades});
    grp.put_attr("num_points_blade", std::vector<int>{meta.num_pts_blade});
    grp.put_attr("num_points_tower", std::vector<int>{meta.num_pts_tower});
    grp.put_attr(
        "rotor_diameter", std::vector<amrex::Real>{info.rotor_diameter});
    grp.put_attr("hub_height", std::vector<amrex::Real>{info.hub_height});
    // clang-format off
    grp.put_attr("base_location", std::vector<amrex::Real>{
            info.base_pos.x(), info.base_pos.y(), info.base_pos.z()});
    // clang-format on

    const size_t nfpts = grid.force.size();
    const size_t nvpts = grid.vel.size();
    grp.def_dim(np_name, nfpts);
    grp.def_dim(nvel_name, nvpts);
    grp.def_var("time", NC_DOUBLE, {nt_name});
    grp.def_var("chord", NC_DOUBLE, {np_name});
    grp.def_var("epsilon", NC_DOUBLE, {np_name, "ndim"});
    grp.def_var("rot_center", NC_DOUBLE, {nt_name, "ndim"});
    grp.def_var("rotor_frame", NC_DOUBLE, {nt_name, "mat_dim"});
    grp.def_var("xyz", NC_DOUBLE, {nt_name, np_name, "ndim"});
    grp.def_var("force", NC_DOUBLE, {nt_name, np_name, "ndim"});
    grp.def_var("orientation", NC_DOUBLE, {nt_name, np_name, "mat_dim"});
    grp.def_var("vel_xyz", NC_DOUBLE, {nt_name, nvel_name, "ndim"});
    grp.def_var("vel", NC_DOUBLE, {nt_name, nvel_name, "ndim"});
    ncf.exit_def_mode();

    {
        auto chord = grp.var("chord");
        chord.put(meta.chord.data(), {0}, {nfpts});
        auto eps = grp.var("epsilon");
        eps.put(grid.epsilon[0].data(), {0, 0}, {nfpts, AMREX_SPACEDIM});
    }

#else
    amrex::ignore_unused(ncfile, meta, info, grid);
#endif
}

void write_netcdf(
    const std::string& ncfile,
    const TurbineBaseData& meta,
    const TurbineInfo& info,
    const ActGrid& grid,
    const amrex::Real time)
{
#ifdef AMR_WIND_USE_NETCDF
    // Only root process handles I/O
    if (!info.is_root_proc) {
        return;
    }

    auto ncf = ncutils::NCFile::open(ncfile, NC_WRITE);
    const std::string nt_name = "num_time_steps";
    // Index of next timestep
    const size_t nt = ncf.dim(nt_name).len();
    const size_t nfpts = grid.force.size();
    const size_t nvpts = grid.vel.size();

    auto grp = ncf.group(info.label);
    grp.var("time").put(&time, {nt}, {1});
    grp.var("rot_center").put(meta.rot_center.data(), {nt, 0}, {1, 3});
    grp.var("rotor_frame").put(meta.rotor_frame.data(), {nt, 0}, {1, 9});
    grp.var("xyz").put(
        grid.pos[0].data(), {nt, 0, 0}, {1, nfpts, AMREX_SPACEDIM});
    grp.var("force").put(
        grid.force[0].data(), {nt, 0, 0}, {1, nfpts, AMREX_SPACEDIM});
    grp.var("orientation")
        .put(grid.orientation[0].data(), {nt, 0, 0}, {1, nfpts, 9});
    grp.var("vel_xyz").put(
        grid.vel_pos[0].data(), {nt, 0, 0}, {1, nvpts, AMREX_SPACEDIM});
    grp.var("vel").put(
        grid.vel[0].data(), {nt, 0, 0}, {1, nvpts, AMREX_SPACEDIM});
#else
    amrex::ignore_unused(ncfile, meta, info, grid, time);
#endif
}

} // namespace amr_wind::actuator::utils
