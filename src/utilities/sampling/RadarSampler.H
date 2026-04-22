#ifndef RADARSAMPLER_H
#define RADARSAMPLER_H

#include <cstdint>
#include "amr-wind/utilities/sampling/SamplerBase.H"
#include "amr-wind/utilities/sampling/SamplingUtils.H"
#include "amr-wind/utilities/sampling/LineSampler.H"
#include "amr-wind/core/vs/vector_space.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::sampling {

/** Sample data along a line
 *  \ingroup sampling
 *
 *  Defines probes along a line defined by `start` and `end` coordinates and
 *  divides it into equal segments as determined by `num_points`. The length of
 *  each segment is given by \f[ (end - start) / (num_points - 1) \f]
 */
class RadarSampler : public SamplerBase::Register<RadarSampler>
{
public:
    static std::string identifier() { return "RadarSampler"; }

    explicit RadarSampler(const CFDSim& /*sim*/);

    ~RadarSampler() override;

    using LosRotType = amrex::Vector<vs::Tensor>;
    using LosUnitType = amrex::Vector<vs::Vector>;

    /** Read user inputs and initialize the sampling object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    enum class phase : std::uint8_t {
        FORWARD,
        FORWARD_PAUSE,
        REVERSE,
        REVERSE_PAUSE
    };

    //! Check and fix the bounds of the sampler so the probes are in the domain
    void check_bounds() override;

    void new_cone();

    [[nodiscard]] amrex::Real total_sweep_time() const;
    amrex::Real periodic_time();
    [[nodiscard]] int sweep_count() const;
    [[nodiscard]] phase determine_operation_phase() const;
    [[nodiscard]] amrex::Real determine_current_sweep_angle() const;

    //! Populate and return a vector of probe locations to be sampled
    void sampling_locations(SampleLocType& /*sample_locs**/) const override;

    //! Populate the vector with coordinates of the sampling locations inside
    //! a box
    void sampling_locations(
        SampleLocType& /*sample_locs*/,
        const amrex::Box& /*box*/) const override;

    bool update_sampling_locations() override;
    void cone_axis_locations(SampleLocType& /*sample_locs*/) const;
    void output_locations(SampleLocType& sample_locs) const override
    {
        cone_axis_locations(sample_locs);
    }

    void post_sample_actions() override;

    void
    define_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void
    populate_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void output_netcdf_data(
        const ncutils::NCGroup& /*unused*/, size_t /*unused*/) const override;

    bool output_netcdf_field(
        const std::vector<amrex::Real>& output_buffer,
        ncutils::NCGroup& grp,
        size_t nt) override;

    //! Name of this sampling object
    [[nodiscard]] std::string label() const override { return m_label; }
    std::string& label() override { return m_label; }

    //! Type of this sampling object
    [[nodiscard]] std::string sampletype() const override
    {
        return identifier();
    }

    //! Unique identifier for this set of probe locations
    [[nodiscard]] int id() const override { return m_id; }
    int& id() override { return m_id; }

    //! Number of sampling container locations in a timestep
    [[nodiscard]] long num_points() const override
    {
        return m_interp_steps * m_ntotal *
               (m_npts + m_ntheta * m_npts * (m_nphi - 1));
    }

    //! Number of probe locations to output
    [[nodiscard]] long num_output_points() const override
    {
        return m_ntotal * m_npts;
    }

    //! Number of sampling container locations per interp step
    [[nodiscard]] long num_points_scan() const
    {
        return m_ntotal * (m_npts + m_ntheta * m_npts * (m_nphi - 1));
    }

    //! Number of probe locations on Radar axis
    [[nodiscard]] long num_points_axis() const { return m_npts; }

    //! Number of probe locations in a cone
    [[nodiscard]] long num_points_cone() const
    {
        return m_npts + (m_ntheta * m_npts * (m_nphi - 1));
    }

    //! Number of probe locations in a spherical cap
    [[nodiscard]] long num_points_quad() const
    {
        return 1 + (m_ntheta * (m_nphi - 1));
    }

    static void line_average(
        const std::vector<amrex::Real>& weights,
        const std::vector<amrex::Real>& values,
        std::vector<amrex::Real>& reduced,
        long offset);

    //! Run data modification for specific sampler
    bool do_data_modification() override { return true; }
    bool do_convert_velocity_los() override { return true; }
    bool do_subsampling_interp() override { return true; }

    //! Modify sample buffer after sampling happens
    std::vector<amrex::Real> modify_sample_data(
        const std::vector<amrex::Real>& sample_data,
        const std::string& /*unused*/) override;

    void calc_lineofsight_velocity(
        const std::vector<std::vector<amrex::Real>>& /*velocity_raw*/,
        int /*interp_idx*/) override;

protected:
    const CFDSim& m_sim;

    std::string m_label;
    int m_id{-1};

    amrex::Real m_radtodeg{180.0_rt / std::numbers::pi_v<amrex::Real>};

    amrex::Vector<amrex::Real> m_start;
    amrex::Vector<amrex::Real> m_end;

    SampleLocType::LocType m_initial_cone;
    SampleLocType::LocType m_current_cones;
    SampleLocType::LocType m_prior_cones;
    SampleLocType::LocType m_sample_cones;

    LosRotType m_los_proj;
    LosUnitType m_los_unit;

    std::vector<amrex::Real> m_los_velocity;
    std::vector<amrex::Real> m_los_velocity_next;
    std::vector<amrex::Real> m_los_velocity_prior;
    std::vector<amrex::Real> m_los_velocity_interp;

    std::vector<amrex::Real> m_weights;
    std::vector<vs::Vector> m_rays;

    amrex::Real m_sample_freq;       // Simulation sample rate
    amrex::Real m_radar_sample_freq; // Actual device sample rate
    int m_npts{0};
    amrex::Real m_cone_angle;
    std::string m_radar_quad_type;
    int m_npts_azimuth{0};
    amrex::Real m_beam_length;
    amrex::Real m_angular_speed;
    amrex::Vector<amrex::Real> m_axis;
    amrex::Vector<amrex::Real> m_vertical;
    amrex::Real m_sweep_angle;
    amrex::Real m_reset_time;
    amrex::Vector<amrex::Real> m_elevation_angles;

    int m_nphi{5};
    int m_ntheta{5};

    vs::Vector m_radar_rotation_axis;

    long m_radar_iter{0};     // Iteration tracking from radar perspective
    amrex::Real m_radar_time; // Current radar time
    amrex::Real
        m_periodic_time; // Running cycle time (Forward, pause, reverse, pause)
    int m_ns;            // Number of radar steps in current timestep
    long m_ntotal;       // Max number of radar timesteps in any timestep
    long m_cone_size;
    phase m_current_phase;

    bool m_debug_print{false};
    bool m_output_cone_points{false};

    long m_interp_steps{2};
    int m_ns_prior;
    long m_ntotal_prior;
    long m_num_points_prior;
    long m_num_output_points_prior;

    amrex::Real m_fill_val{-99999.99_rt};
};

} // namespace amr_wind::sampling

#endif /* RADARSAMPLER_H */
