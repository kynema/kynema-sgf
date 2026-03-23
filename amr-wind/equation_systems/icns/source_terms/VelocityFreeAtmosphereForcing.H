#ifndef VELOCITYFREEATMOSPHEREFORCING_H
#define VELOCITYFREEATMOSPHEREFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::pde::icns {

class VelocityFreeAtmosphereForcing
    : public MomentumSource::Register<VelocityFreeAtmosphereForcing>
{
public:
    static std::string identifier() { return "VelocityFreeAtmosphereForcing"; }

    explicit VelocityFreeAtmosphereForcing(const CFDSim& sim);

    ~VelocityFreeAtmosphereForcing() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState /*fstate*/,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;
    std::string m_1d_rans_filename;
    const Field& m_velocity;
    amrex::Vector<amrex::Real> m_wind_heights;
    amrex::Vector<amrex::Real> m_u_values;
    amrex::Vector<amrex::Real> m_v_values;
    amrex::Vector<amrex::Real> m_w_values;
    amrex::Gpu::DeviceVector<amrex::Real> m_wind_heights_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_u_values_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_v_values_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_w_values_d;
    amrex::Real m_meso_start{1e30};
    amrex::Real m_meso_timescale{30};
    const CFDSim& m_sim;
};

} // namespace amr_wind::pde::icns
#endif
