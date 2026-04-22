#ifndef PERTURBATIONFORCING_H
#define PERTURBATIONFORCING_H

#include "amr-wind/equation_systems/temperature/TemperatureSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/transport_models/TransportModel.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::temperature {

class PerturbationForcing
    : public TemperatureSource::Register<PerturbationForcing>
{
public:
    static std::string identifier() { return "PerturbationForcing"; }

    explicit PerturbationForcing(const CFDSim& sim);

    ~PerturbationForcing() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState /*fstate*/,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const CFDSim& m_sim;
    const SimTime& m_time;
    const amrex::AmrCore& m_mesh;
    amrex::Vector<amrex::Real> m_start;
    amrex::Vector<amrex::Real> m_end;
    int m_time_index{100};
    int m_start_level{0};
    amrex::Real m_pert_amplitude{1.0_rt};
};

} // namespace amr_wind::pde::temperature
#endif
