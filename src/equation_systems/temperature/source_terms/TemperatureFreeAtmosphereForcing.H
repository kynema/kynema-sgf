#ifndef TEMPERATUREFREEATMOSPHEREFORCING_H
#define TEMPERATUREFREEATMOSPHEREFORCING_H

#include "src/equation_systems/temperature/TemperatureSource.H"
#include "src/core/SimTime.H"
#include "src/CFDSim.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::temperature {

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
    amrex::Real m_meso_start{1e30_rt};
    amrex::Real m_meso_timescale{30.0_rt};
    bool m_horizontal_sponge{false};
    amrex::Real m_sponge_strength{1.0_rt};
    amrex::Real m_sponge_density{1.0_rt};
    amrex::Real m_sponge_distance_west{0.0_rt};
    amrex::Real m_sponge_distance_east{0.0_rt};
    amrex::Real m_sponge_distance_south{0.0_rt};
    amrex::Real m_sponge_distance_north{0.0_rt};
    bool m_sponge_west{false};
    bool m_sponge_east{false};
    bool m_sponge_south{false};
    bool m_sponge_north{false};
    const CFDSim& m_sim;
};

} // namespace kynema_sgf::pde::temperature
#endif
