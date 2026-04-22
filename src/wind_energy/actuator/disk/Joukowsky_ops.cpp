#include "amr-wind/wind_energy/actuator/disk/Joukowsky_ops.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/wind_energy/actuator/disk/disk_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator::ops::joukowsky {
void check_for_parse_conflicts(const utils::ActParser& pp)
{
    auto error_collector = ops::base::check_for_parse_conflicts(pp);
    // clang-format off
    base::collect_parse_dependencies(pp, "thrust_coeff", "rpm", error_collector);
    // clang-format on
    ops::base::check_error_stream(error_collector);
}

void optional_parameters(JoukowskyData& meta, const utils::ActParser& pp)
{
    ops::base::optional_parameters(meta, pp);
    pp.query("use_root_correction", meta.use_root_correction);
    pp.query("use_tip_correction", meta.use_tip_correction);
    if (pp.contains("vortex_core_size")) {
        pp.get("vortex_core_size", meta.vortex_core_size);
    } else {
        // Default to 20% of radius
        meta.vortex_core_size = 0.1_rt * meta.diameter;
    }
    pp.query("root_correction_exponent", meta.root_correction_exponent);
    pp.query("root_correction_coefficient", meta.root_correction_coefficient);
    pp.query("num_blades", meta.num_blades);
    pp.query("ct_region2", meta.Ct_rated);
    pp.query("S0_alpha1", meta.S0_alpha1);
    pp.query("S0_alpha2", meta.S0_alpha2);
    pp.query("awc_amplitude", meta.awc_amplitude);
    pp.query("awc_azimuthal_mode", meta.awc_azimuthal_mode);
    pp.query("awc_angular_frequency", meta.awc_angular_frequency);
    pp.query("awc_clocking_angle", meta.awc_clocking_angle);
}

void required_parameters(JoukowskyData& meta, const utils::ActParser& pp)
{
    ops::base::required_parameters(meta, pp);
    pp.get("num_points_t", meta.num_vel_pts_t);
    pp.get("num_points_r", meta.num_vel_pts_r);
    pp.getarr("rpm", meta.angular_velocity);
    // Convert from rpm to rad/s
    for (amrex::Real& i : meta.angular_velocity) {
        i *= std::numbers::pi_v<amrex::Real> / 30.0_rt;
    }
}

void parse_and_gather_params(const utils::ActParser& pp, JoukowskyData& data)
{
    check_for_parse_conflicts(pp);
    required_parameters(data, pp);
    optional_parameters(data, pp);
    ops::base::final_checks(data);
    data.num_force_pts = data.num_vel_pts_r * data.num_vel_pts_t;
    data.num_vel_pts = data.num_force_pts * 2;
    data.dr = 0.5_rt * data.diameter / data.num_vel_pts_r;
}

void update_disk_points(Joukowsky::DataType& data)
{
    auto& grid = data.grid();
    auto& meta = data.meta();

    base::compute_and_normalize_coplanar_vector(meta);
    data.info().bound_box = base::compute_bounding_box(meta);

    const auto& sVec = meta.sample_vec;
    const auto& nVec = meta.normal_vec;

    // force points
    base::compute_disk_points(meta, grid.pos, nVec, 0, 0);
    // velocity points upstream
    base::compute_disk_points(
        meta, grid.vel_pos, sVec, 0, meta.diameters_to_sample);
    // velocity points at the disk
    base::compute_disk_points(
        meta, grid.vel_pos, nVec, meta.num_vel_pts / 2, 0);
}

void prepare_netcdf_file(
    const std::string& name,
    const JoukowskyData& data,
    const ActInfo& info,
    const ActGrid& grid)
{
#ifdef AMR_WIND_USE_NETCDF
    using dvec = std::vector<amrex::Real>;
    // Only root process handles I/O
    if (info.root_proc != amrex::ParallelDescriptor::MyProc()) {
        return;
    }
    auto ncf = ncutils::NCFile::create(name, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    const std::string np_name = "num_actuator_points";
    const std::string nv_name = "num_velocity_points";
    const std::string nr_name = "num_radial_points";

    ncf.enter_def_mode();
    ncf.put_attr("title", "AMR-Wind ActuatorDisk actuator output [Joukowsky]");
    ncf.put_attr("version", ioutils::amr_wind_version());
    ncf.put_attr("created_on", ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim("ndim", AMREX_SPACEDIM);
    ncf.def_dim(nr_name, data.num_vel_pts_r);

    auto grp = ncf.def_group(info.label);
    // clang-format off
    grp.put_attr("normal",
        dvec{data.normal_vec.x(), data.normal_vec.y(), data.normal_vec.z()});
    grp.put_attr("sample_normal",
        dvec{data.sample_vec.x(), data.sample_vec.y(), data.sample_vec.z()});
    // clang-format on
    grp.put_attr("diameter", dvec{data.diameter});
    grp.put_attr("epsilon", dvec{data.epsilon});
    grp.put_attr("sample_diameters", dvec{data.diameters_to_sample});
    grp.def_dim(np_name, data.num_force_pts);
    grp.def_dim(nv_name, data.num_vel_pts);
    grp.def_var("time", NC_DOUBLE, {nt_name});
    grp.def_var("xyz", NC_DOUBLE, {np_name, "ndim"});
    grp.def_var("xyz_v", NC_DOUBLE, {nv_name, "ndim"});
    grp.def_var("vref", NC_DOUBLE, {nt_name, "ndim"});
    grp.def_var("vdisk", NC_DOUBLE, {nt_name, "ndim"});
    grp.def_var("tsr", NC_DOUBLE, {nt_name});
    grp.def_var("ct", NC_DOUBLE, {nt_name});
    grp.def_var("cp", NC_DOUBLE, {nt_name});
    grp.def_var("power", NC_DOUBLE, {nt_name});
    grp.def_var("density", NC_DOUBLE, {nt_name});
    grp.def_var("total_disk_force", NC_DOUBLE, {nt_name, "ndim"});
    grp.def_var("angular_velocity", NC_DOUBLE, {nt_name});
    grp.def_var("f_normal", NC_DOUBLE, {nt_name, nr_name});
    grp.def_var("f_theta", NC_DOUBLE, {nt_name, np_name});
    ncf.exit_def_mode();

    {
        {
            const auto npts = static_cast<size_t>(data.num_force_pts);
            const std::vector<size_t> start{0, 0};
            const std::vector<size_t> count{npts, AMREX_SPACEDIM};
            grp.var("xyz").put(grid.pos[0].data(), start, count);
        }
        {
            const auto npts = static_cast<size_t>(data.num_vel_pts);
            const std::vector<size_t> start{0, 0};
            const std::vector<size_t> count{npts, AMREX_SPACEDIM};
            grp.var("xyz_v").put(grid.vel_pos[0].data(), start, count);
        }
    }
#else
    amrex::ignore_unused(name, data, info, grid);
#endif
}

void write_netcdf(
    const std::string& name,
    const JoukowskyData& data,
    const ActInfo& info,
    const ActGrid& /*unused*/,
    const amrex::Real time)
{
#ifdef AMR_WIND_USE_NETCDF
    // Only root process handles I/O
    if (info.root_proc != amrex::ParallelDescriptor::MyProc()) {
        return;
    }
    auto ncf = ncutils::NCFile::open(name, NC_WRITE);
    const std::string nt_name = "num_time_steps";
    // Index of next timestep
    const size_t nt = ncf.dim(nt_name).len();
    const size_t nr = data.num_vel_pts_r;
    const size_t np = data.num_force_pts;
    auto grp = ncf.group(info.label);
    grp.var("time").put(&time, {nt}, {1});
    grp.var("vref").put(
        data.reference_velocity.data(), {nt, 0}, {1, AMREX_SPACEDIM});
    grp.var("vdisk").put(
        data.mean_disk_velocity.data(), {nt, 0}, {1, AMREX_SPACEDIM});
    grp.var("tsr").put(&data.current_tip_speed_ratio, {nt}, {1});
    grp.var("ct").put(&data.current_ct, {nt}, {1});
    grp.var("cp").put(&data.current_cp, {nt}, {1});
    grp.var("power").put(&data.current_power, {nt}, {1});
    grp.var("density").put(&data.density, {nt}, {1});
    grp.var("total_disk_force")
        .put(data.disk_force.data(), {nt, 0}, {1, AMREX_SPACEDIM});
    grp.var("angular_velocity").put(&data.current_angular_velocity, {nt}, {1});
    grp.var("f_normal").put(data.f_normal.data(), {nt, 0}, {1, nr});
    grp.var("f_theta").put(data.f_theta.data(), {nt, 0}, {1, np});
#else
    amrex::ignore_unused(name, data, info, time);
#endif
}

} // namespace amr_wind::actuator::ops::joukowsky
