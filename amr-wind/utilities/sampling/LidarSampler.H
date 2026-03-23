#ifndef LIDARSAMPLER_H
#define LIDARSAMPLER_H

#include "amr-wind/utilities/sampling/LineSampler.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind::sampling {

/** Sample data along a line that moves in time
 *  \ingroup sampling
 *
 *  Defines probes along a line defined by a time history from `m_azimuth_table`
 *  and `m_elevation_table` with a length of `m_length`. This determines the
 * `start` and `end` points at every time-step and divides the line into equal
 * segments as determined by `num_points`. The length of each segment is given
 * by \f[ m_length / (num_points - 1) \f]
 */
class LidarSampler : public LineSampler
{
public:
    // cppcheck-suppress duplInheritedMember
    static std::string identifier() { return "LidarSampler"; }

    explicit LidarSampler(const CFDSim& /*sim*/);

    /** Read user inputs and initialize the sampling object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    /** Update the position of the sampling locations in time
     *
     *
     */
    bool update_sampling_locations() override;

    void post_sample_actions() override {};

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

protected:
    amrex::Vector<amrex::Real> m_origin;
    amrex::Vector<amrex::Real> m_time_table;
    amrex::Vector<amrex::Real> m_azimuth_table;
    amrex::Vector<amrex::Real> m_elevation_table;
    amrex::Real m_length{0};
    bool m_periodic{true};
    amrex::Real m_period{std::numeric_limits<amrex::Real>::max()};
};

} // namespace amr_wind::sampling

#endif /* LIDARSAMPLER_H */
