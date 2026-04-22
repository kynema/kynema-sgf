#ifndef PROBESAMPLER_H
#define PROBESAMPLER_H

#include "amr-wind/utilities/sampling/SamplerBase.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::sampling {

/** Sample data on a cloud of points
 *  \ingroup sampling
 *
 *  Reads a cloud of point locations from a text file and samples data on those
 *  points.
 */
class ProbeSampler : public SamplerBase::Register<ProbeSampler>
{
public:
    static std::string identifier() { return "ProbeSampler"; }

    //! Type of this sampling object
    [[nodiscard]] std::string sampletype() const override
    {
        return identifier();
    }

    explicit ProbeSampler(const CFDSim& /*unused*/);

    ~ProbeSampler() override;

    void initialize(const std::string& key) override;

    //! Check and fix the bounds of the sampler so the probes are in the domain
    void check_bounds() override;

    //! Populate the vector with coordinates of the sampling locations
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

    [[nodiscard]] std::string label() const override { return m_label; }
    std::string& label() override { return m_label; }
    [[nodiscard]] int id() const override { return m_id; }
    int& id() override { return m_id; }

    [[nodiscard]] long num_points() const override { return m_npts; }

    //! Number of output probe locations (after data reduction etc.)
    [[nodiscard]] long num_output_points() const override { return m_npts; }

private:
    const CFDSim& m_sim;
    SampleLocType m_probes;

    amrex::Vector<amrex::Real> m_offset_vector{0.0_rt, 0.0_rt, 0.0_rt};
    amrex::Vector<amrex::Real> m_poffsets;

    std::string m_label;
    int m_id{-1};
    int m_npts{0};
};

} // namespace amr_wind::sampling

#endif /* PROBESAMPLER_H */
