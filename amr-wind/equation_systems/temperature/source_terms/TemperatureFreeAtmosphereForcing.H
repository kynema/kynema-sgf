#ifndef TEMPERATUREFREEATMOSPHEREFORCING_H
#define TEMPERATUREFREEATMOSPHEREFORCING_H

#include "amr-wind/equation_systems/temperature/TemperatureSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::temperature {

class TemperatureFreeAtmosphereForcing
    : public TemperatureSource::Register<TemperatureFreeAtmosphereForcing>
{
public:
    static std::string identifier()
    {
        return "TemperatureFreeAtmosphereForcing";
    }

    explicit TemperatureFreeAtmosphereForcing(const CFDSim& sim);

    ~TemperatureFreeAtmosphereForcing() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState /*fstate*/,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const amrex::AmrCore& m_mesh;
    const Field& m_temperature;
    amrex::Vector<amrex::Real> m_theta_heights;
    amrex::Vector<amrex::Real> m_theta_values;
    amrex::Gpu::DeviceVector<amrex::Real> m_theta_heights_d;
    amrex::Gpu::DeviceVector<amrex::Real> m_theta_values_d;
    amrex::Real m_meso_start{600.0_rt};
    amrex::Real m_meso_timescale{30.0_rt};
    bool m_horizontal_sponge{false};
    amrex::Real m_sponge_strength{1.0_rt};
    amrex::Real m_sponge_density{1.0_rt};
    amrex::Real m_sponge_distance_west{-1000.0_rt};
    amrex::Real m_sponge_distance_east{1000.0_rt};
    amrex::Real m_sponge_distance_south{-1000.0_rt};
    amrex::Real m_sponge_distance_north{1000.0_rt};
    int m_sponge_west{0};
    int m_sponge_east{1};
    int m_sponge_south{0};
    int m_sponge_north{1};
    const CFDSim& m_sim;
};

} // namespace amr_wind::pde::temperature
#endif
