#ifndef DTUSPINNERSAMPLER_H
#define DTUSPINNERSAMPLER_H

#include "amr-wind/utilities/sampling/SamplingUtils.H"
#include "amr-wind/utilities/sampling/LidarSampler.H"
#include "amr-wind/core/vs/vector_space.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::sampling {

struct PrismParameters
{
    amrex::Real theta0{90.0_rt};
    amrex::Real rot{3.5_rt};
    amrex::Real azimuth{15.2_rt};
};

/** Sample data along a line that moves in time
 *  \ingroup sampling
 *
 *  Defines probes along a line defined by a time history from `m_azimuth_table`
 *  and `m_elevation_table` with a length of `m_length`. This determines the
 * `start` and `end` points at every time-step and divides the line into equal
 * segments as determined by `num_points`. The length of each segment is given
 * by \f[ m_length / (num_points - 1) \f]
 */
class DTUSpinnerSampler : public LidarSampler
{
public:
    // cppcheck-suppress duplInheritedMember
    static std::string identifier() { return "DTUSpinnerSampler"; }

    explicit DTUSpinnerSampler(const CFDSim& /*sim*/);

    /** Read user inputs and initialize the sampling object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    /** Update the position of the sampling locations in time
     *
     *
     */
    void sampling_locations(SampleLocType& /*sample_locs*/) const override;

    //! Populate the vector with coordinates of the sampling locations inside
    //! a box
    void sampling_locations(
        SampleLocType& /*sample_locs*/,
        const amrex::Box& /*box*/) const override;

    static vs::Vector generate_lidar_pattern(
        PrismParameters InnerPrism,
        PrismParameters OuterPrism,
        amrex::Real time);

    static vs::Vector adjust_lidar_pattern(
        vs::Vector beamPt,
        amrex::Real yaw,
        amrex::Real pitch,
        amrex::Real roll);

#ifdef AMR_WIND_USE_OPENFAST

    void bcast_turbine(
        amrex::Array<amrex::Real, 18>& /*turbine_pack*/, int /*root_proc*/);
    void get_turbine_data(const std::string& turbine_label);

#endif

    /** Determine how the subsampling will be done
     *
     */
    bool update_sampling_locations() override;

    void post_sample_actions() override {};

    [[nodiscard]] long num_points() const override
    {
        return (m_beam_points * m_ntotal);
    }

    [[nodiscard]] long num_output_points() const override
    {
        return (m_beam_points * m_ntotal);
    }

    //! Type of this sampling object
    [[nodiscard]] std::string sampletype() const override
    {
        return identifier();
    }

    void
    define_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void
    populate_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void output_netcdf_data(
        const ncutils::NCGroup& /*unused*/, size_t /*unused*/) const override;

private:
    // Number of subsamples in this timestep
    int m_ns{1};

    // Max number of subsamples in all timesteps
    long m_ntotal{1};

    int m_update_count{0};
    const amrex::Real m_radtodeg{180.0_rt / std::numbers::pi_v<amrex::Real>};
    const amrex::Real m_pi{std::numbers::pi_v<amrex::Real>};
    const amrex::Real m_twopi{2.0_rt * std::numbers::pi_v<amrex::Real>};
    const amrex::Real m_threepi{3.0_rt * std::numbers::pi_v<amrex::Real>};
    amrex::Real m_time_sampling{0.0_rt};
    PrismParameters m_InnerPrism;
    PrismParameters m_OuterPrism;
    amrex::Vector<amrex::Real> m_lidar_center{0.0_rt, 0.0_rt, 0.0_rt};
    amrex::Real m_scan_time{2.0_rt};
    amrex::Real m_num_samples{984.0_rt};
    amrex::Real m_beam_length{270.0_rt};
    int m_beam_points{432};
    amrex::Real m_fixed_yaw{0.0_rt};
    amrex::Real m_fixed_roll{0.0_rt};
    amrex::Real m_fixed_tilt{0.0_rt};
    amrex::Real m_hub_yaw{0.0_rt};
    amrex::Real m_hub_roll{0.0_rt};
    amrex::Real m_hub_tilt{0.0_rt};
    std::string m_spinner_mode{"fixed"};
    std::string m_turbine_label;
    bool m_hub_debug{false};
    vs::Vector m_hub_location{0, 0, 0};
    vs::Vector m_hub_location_init{0, 0, 0};

#ifdef AMR_WIND_USE_OPENFAST

    amrex::Array<float, 3> m_turbine_base_pos = {0, 0, 0};
    amrex::Array<float, 3> m_current_hub_abs_pos = {0, 0, 0};
    amrex::Array<float, 3> m_current_hub_rot_vel = {0, 0, 0};
    amrex::Array<amrex::Real, 9> m_current_hub_orient = {0, 0, 0, 0, 0,
                                                         0, 0, 0, 0};

#endif

    amrex::Real m_last_hub_yaw{0.0_rt};
    amrex::Real m_last_hub_roll{0.0_rt};
    amrex::Real m_last_hub_tilt{0.0_rt};
    amrex::Vector<amrex::Real> m_last_lidar_center{0.0_rt, 0.0_rt, 0.0_rt};
    amrex::Vector<amrex::Real> m_step_lidar_center{0.0_rt, 0.0_rt, 0.0_rt};

    amrex::Real m_fill_val{-99999.99_rt};
};

} // namespace amr_wind::sampling

#endif /* DTUSPINNERSAMPLER_H */
