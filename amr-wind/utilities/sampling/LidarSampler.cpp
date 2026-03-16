#include "amr-wind/utilities/sampling/LidarSampler.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {
namespace sampling {

LidarSampler::LidarSampler(const CFDSim& sim) : LineSampler(sim) {}

void LidarSampler::initialize(const std::string& key)
{
    // Read in new inputs specific to this class
    amrex::ParmParse pp(key);

    // This is the origin of the scan (x, y, z) [m]
    pp.getarr("origin", m_origin);
    AMREX_ALWAYS_ASSERT(static_cast<int>(m_origin.size()) == AMREX_SPACEDIM);

    // The number of points
    pp.get("num_points", m_npts);

    // The length of the beam [m]
    pp.get("length", m_length);

    // The time table [s]
    pp.getarr("time_table", m_time_table);
    // Azimuth angle (this is in the N, E, S, W direction) [degrees]
    pp.getarr("azimuth_table", m_azimuth_table);
    // Elevation angle [degrees]
    pp.getarr("elevation_table", m_elevation_table);

    // Establish if the time table is periodic (repeats in time) or not
    pp.query("periodic", m_periodic);

    // Number of elements in the table
    int np = static_cast<int>(m_time_table.size());

    // Ensure that the tables have the same size
    if (m_azimuth_table.size() != np) {
        amrex::Abort(
            "azimuth_table must have same number of entries as time_table ");
    }
    if (m_elevation_table.size() != np) {
        amrex::Abort(
            "elevation_table must have same number of entries as time_table ");
    }

    // Initialize the star and end variables
    m_start = m_origin;
    m_end = {0.0_rt, 0.0_rt, 0.0_rt};

    // The period to know if the table repeats in time
    m_period = m_time_table[np - 1] - m_time_table[0];

    LidarSampler::update_sampling_locations();

    check_bounds();
}

bool LidarSampler::update_sampling_locations()
{

    amrex::Real time = m_sim.time().current_time();

    // If the time table is periodic, make sure to get the correct time
    if (m_periodic) {
        time = std::fmod(time, m_period);
    }

    // The current azimuth angle
    const amrex::Real current_azimuth = utils::radians(
        ::amr_wind::interp::linear(m_time_table, m_azimuth_table, time));

    const amrex::Real current_elevation = utils::radians(
        90.0_rt -
        ::amr_wind::interp::linear(m_time_table, m_elevation_table, time));

    // Coordinate transform spherical to cartesian
    m_end[0] = m_start[0] + (m_length * std::cos(current_azimuth) *
                             std::sin(current_elevation));
    m_end[1] = m_start[1] + (m_length * std::sin(current_azimuth) *
                             std::sin(current_elevation));
    m_end[2] = m_start[2] + (m_length * std::cos(current_elevation));

    return true;
}

#ifdef AMR_WIND_USE_NETCDF
void LidarSampler::define_netcdf_metadata(const ncutils::NCGroup& grp) const
{
    grp.put_attr("sampling_type", identifier());
    grp.put_attr("start", m_start);
    grp.put_attr("end", m_end);
    grp.put_attr("time_table", m_time_table);
    grp.put_attr("azimuth_table", m_azimuth_table);
    grp.put_attr("elevation_table", m_elevation_table);
    grp.def_var("points", NC_DOUBLE, {"num_time_steps", "num_points", "ndim"});
}

void LidarSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
void LidarSampler::output_netcdf_data(
    const ncutils::NCGroup& grp, const size_t nt) const
{
    // Write the coordinates every time
    std::vector<size_t> start{nt, 0, 0};
    std::vector<size_t> count{1, 0, AMREX_SPACEDIM};
    SampleLocType sample_locs;
    sampling_locations(sample_locs);
    auto xyz = grp.var("points");
    count[1] = num_points();
    const auto& locs = sample_locs.locations();
    xyz.put(locs[0].begin(), start, count);
}
#else
void LidarSampler::define_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
void LidarSampler::populate_netcdf_metadata(
    const ncutils::NCGroup& /*unused*/) const
{}
void LidarSampler::output_netcdf_data(
    const ncutils::NCGroup& /*unused*/, const size_t /*unused*/) const
{}
#endif

} // namespace sampling

template struct ::amr_wind::sampling::SamplerBase::Register<
    ::amr_wind::sampling::LidarSampler>;

} // namespace amr_wind
