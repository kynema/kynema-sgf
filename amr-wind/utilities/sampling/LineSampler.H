#ifndef LINESAMPLER_H
#define LINESAMPLER_H

#include "amr-wind/utilities/sampling/SamplerBase.H"

namespace amr_wind::sampling {

/** Sample data along a line
 *  \ingroup sampling
 *
 *  Defines probes along a line defined by `start` and `end` coordinates and
 *  divides it into equal segments as determined by `num_points`. The length of
 *  each segment is given by \f[ (end - start) / (num_points - 1) \f]
 */
class LineSampler : public SamplerBase::Register<LineSampler>
{
public:
    static std::string identifier() { return "LineSampler"; }

    explicit LineSampler(const CFDSim& /*sim*/);

    ~LineSampler() override;

    /** Read user inputs and initialize the sampling object
     *
     *  \param key Prefix used to parse inputs from file
     */
    void initialize(const std::string& key) override;

    //! Check and fix the bounds of the sampler so the probes are in the domain
    void check_bounds() override;

    //! Populate and return a vector of probe locations to be sampled
    void sampling_locations(SampleLocType& /*sample_locs*/) const override;

    //! Populate the vector with coordinates of the sampling locations inside
    //! a box
    void sampling_locations(
        SampleLocType& /*sample_locs*/,
        const amrex::Box& /*box*/) const override;

    void output_locations(SampleLocType& sample_locs) const override
    {
        sampling_locations(sample_locs);
    }

    void
    define_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;
    void
    populate_netcdf_metadata(const ncutils::NCGroup& /*unused*/) const override;

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

    //! Number of probe locations along the line
    [[nodiscard]] long num_points() const override { return m_npts; }

    //! Number of output probe locations (after data reduction etc.)
    [[nodiscard]] long num_output_points() const override { return m_npts; }

protected:
    const CFDSim& m_sim;

    amrex::Vector<amrex::Real> m_start;
    amrex::Vector<amrex::Real> m_end;

    std::string m_label;
    int m_id{-1};
    int m_npts{0};
};

} // namespace amr_wind::sampling

#endif /* LINESAMPLER_H */
