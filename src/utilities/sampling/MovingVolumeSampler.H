#ifndef MOVINGVOLUMESAMPLER_H
#define MOVINGVOLUMESAMPLER_H

#include "src/utilities/sampling/VolumeSampler.H"

namespace kynema_sgf::actuator {
class ActuatorModel;
}

namespace kynema_sgf::sampling {

/** Volume sampler expressed in and moving with an actuator reference frame. */
class MovingVolumeSampler : public VolumeSampler
{
public:
    static std::string identifier() { return "MovingVolumeSampler"; }

    explicit MovingVolumeSampler(const CFDSim& sim);

    void initialize(const std::string& key) override;
    void check_bounds() override;
    void sampling_locations(SampleLocType& sample_locs) const override;
    void sampling_locations(
        SampleLocType& sample_locs, const amrex::Box& box) const override;
    bool update_sampling_locations() override;

    [[nodiscard]] std::string sampletype() const override
    {
        return identifier();
    }

    void define_netcdf_metadata(const ncutils::NCGroup& grp) const override;
    void
    output_netcdf_data(const ncutils::NCGroup& grp, size_t nt) const override;
    void populate_info_file(std::ostream& fh) const override;

private:
    struct VolumeGeometry
    {
        amrex::Vector<amrex::Real> origin;
        amrex::Vector<amrex::Real> axis1;
        amrex::Vector<amrex::Real> axis2;
        amrex::Vector<amrex::Real> axis3;
    };

    bool update_geometry(amrex::Real time);

    VolumeGeometry m_reference;
    VolumeGeometry m_geometry;
    const ::kynema_sgf::actuator::ActuatorModel* m_actuator{nullptr};
    std::string m_actuator_label;
    bool m_defer_bounds_check{true};
};

} // namespace kynema_sgf::sampling

#endif
