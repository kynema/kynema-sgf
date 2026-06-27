#ifndef EBDRAGTEMPFORCING_H
#define EBDRAGTEMPFORCING_H

#include "src/equation_systems/temperature/TemperatureSource.H"
#include "src/core/SimTime.H"
#include "src/CFDSim.H"
#include "src/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::temperature {

class EBDragTempForcing : public TemperatureSource::Register<EBDragTempForcing>
{
public:
    static std::string identifier() { return "EBDragTempForcing"; }

    explicit EBDragTempForcing(const CFDSim& sim);

    ~EBDragTempForcing() override;

    void operator()(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override;

private:
    const SimTime& m_time;
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;
    const Field& m_temperature;
    amrex::Real m_drag_coefficient{1.0_rt};

    std::string m_wall_het_model{"none"};
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
    amrex::Real m_kappa{0.41_rt};
    amrex::Real m_gamma_m{5.0_rt};
    amrex::Real m_beta_m{16.0_rt};
    amrex::Real m_gamma_h{5.0_rt};
    amrex::Real m_beta_h{16.0_rt};
    amrex::Real m_soil_temperature{300};
    amrex::Real m_forcing_time_factor{5.0_rt};
    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};

    // IDW soil temperature model parameters
    bool m_use_idw_model{false};
    amrex::Gpu::DeviceVector<amrex::Real> m_idw_x;
    amrex::Gpu::DeviceVector<amrex::Real> m_idw_y;
    amrex::Gpu::DeviceVector<amrex::Real> m_idw_z;
    amrex::Gpu::DeviceVector<amrex::Real> m_idw_temp;
    int m_idw_num_points{0};

    // Heat flux based soil temperature model parameters
    bool m_use_heatflux_model{false};
    amrex::Real m_surface_heatflux{0.0_rt}; // W/m^2
    amrex::Real m_thermal_conductivity{
        0.026_rt}; // W/(m·K) - typical for air at surface
    amrex::Real m_boundary_layer_height{1.0_rt}; // meters
};

} // namespace kynema_sgf::pde::temperature
#endif
