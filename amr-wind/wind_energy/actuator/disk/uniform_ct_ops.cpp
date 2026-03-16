#include "amr-wind/wind_energy/actuator/disk/uniform_ct_ops.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "AMReX_REAL.H"
#include <ostream>

using namespace amrex::literals;

namespace amr_wind::actuator::ops::uniformct {

void check_for_removed_syntax(
    const utils::ActParser& pp, std::ostringstream& stream)
{
    const std::string error_start = "ERROR ActuatorDisk:: ";
    if (pp.contains("num_vel_points_r")) {
        stream << error_start
               << "'num_vel_points_r' has been replaced with 'num_points_r'\n";
    }
    if (pp.contains("num_vel_points_t")) {
        stream << error_start
               << "'num_vel_points_t' has been replaced with 'num_points_t'\n";
    }
    if (pp.contains("num_force_points")) {
        stream << error_start
               << "'num_force_points' has been replaced with 'num_points_r'\n";
    }
    if (pp.contains("num_theta_force_points")) {
        stream << error_start
               << "'num_theta_force_points' has been replaced with "
                  "'num_points_t'\n";
    }
}

void check_for_parse_conflicts(const utils::ActParser& pp)
{
    auto error_collector = ops::base::check_for_parse_conflicts(pp);
    check_for_removed_syntax(pp, error_collector);
    ops::base::collect_parse_dependencies_one_way(
        pp, "num_points_t", "spreading", error_collector);
    ops::base::check_error_stream(error_collector);
}
void optional_parameters(DiskBaseData& meta, const utils::ActParser& pp)
{
    ops::base::optional_parameters(meta, pp);
    // params for the sampling disk
    pp.query("num_points_t", meta.num_vel_pts_t);
    pp.query("spreading_type", meta.spreading_type);
}
void required_parameters(DiskBaseData& meta, const utils::ActParser& pp)
{
    ops::base::required_parameters(meta, pp);
    pp.get("num_points_r", meta.num_force_pts);
}
void parse_and_gather_params(const utils::ActParser& pp, UniformCtData& data)
{
    check_for_parse_conflicts(pp);
    required_parameters(data, pp);
    optional_parameters(data, pp);
    ops::base::final_checks(data);
    data.num_vel_pts_r = data.num_force_pts;
    data.dr = 0.5_rt * data.diameter / data.num_force_pts;
    // 2x 1 for sampling up stream and one for sampling at the disk
    data.num_vel_pts = data.num_vel_pts_r * data.num_vel_pts_t * 2;
}

void update_disk_points(UniformCt::DataType& data)
{
    auto& grid = data.grid();
    auto& meta = data.meta();

    base::compute_and_normalize_coplanar_vector(meta);
    data.info().bound_box = base::compute_bounding_box(meta);

    const auto& cVec = meta.coplanar_vec;
    const auto& sVec = meta.sample_vec;
    const auto& nVec = meta.normal_vec;

    const auto& cc = meta.center;
    {
        const auto& dr = meta.dr;
        for (int i = 0; i < meta.num_force_pts; ++i) {
            grid.pos[i] = cc + (i + 0.5_rt) * dr * cVec;
        }
    }
    base::compute_disk_points(
        meta, grid.vel_pos, sVec, 0, meta.diameters_to_sample);
    base::compute_disk_points(
        meta, grid.vel_pos, nVec, meta.num_vel_pts / 2, 0);
}
} // namespace amr_wind::actuator::ops::uniformct
