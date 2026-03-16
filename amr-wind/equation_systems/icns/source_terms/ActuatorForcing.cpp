#include "amr-wind/equation_systems/icns/source_terms/ActuatorForcing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/wind_energy/actuator/Actuator.H"

#include "AMReX_Gpu.H"

namespace amr_wind::pde::icns {

ActuatorForcing::ActuatorForcing(const CFDSim& sim)
    : m_act_src(sim.repo().get_field("actuator_src_term"))
{
    if (!sim.physics_manager().contains("Actuator")) {
        amrex::Abort("ActuatorForcing requires Actuator physics to be active");
    }
}

ActuatorForcing::~ActuatorForcing() = default;

void ActuatorForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto varr = m_act_src(lev).const_array(mfi);

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        src_term(i, j, k, 0) += varr(i, j, k, 0);
        src_term(i, j, k, 1) += varr(i, j, k, 1);
        src_term(i, j, k, 2) += varr(i, j, k, 2);
    });
}

} // namespace amr_wind::pde::icns
