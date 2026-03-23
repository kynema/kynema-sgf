#ifndef DRAGTEMPFORCING_H
#define DRAGTEMPFORCING_H

#include "amr-wind/equation_systems/temperature/TemperatureSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::temperature {

class DragTempForcing : public TemperatureSource::Register<DragTempForcing>
{
public:
    static std::string identifier() { return "DragTempForcing"; }

    explicit DragTempForcing(const CFDSim& sim);

    ~DragTempForcing() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState /*fstate*/,
        const amrex::Array4<amrex::Real>& src_term) const override;

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
    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};
};

} // namespace amr_wind::pde::temperature
#endif
