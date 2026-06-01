#include "src/equation_systems/icns/source_terms/ActuatorForcing.H"
#include "src/CFDSim.H"
#include "src/wind_energy/actuator/Actuator.H"

#include "AMReX_Gpu.H"

namespace kynema_sgf::pde::icns {

ActuatorForcing::ActuatorForcing(const CFDSim& sim)
    : m_act_src(sim.repo().get_field("actuator_src_term"))
{
    if (!sim.physics_manager().contains("Actuator")) {
        amrex::Abort("ActuatorForcing requires Actuator physics to be active");
    }
}

ActuatorForcing::~ActuatorForcing() = default;

void ActuatorForcing::operator()(
    const int lev, const FieldState /*fstate*/, amrex::MultiFab& src_term) const
{
    auto const& src_arrs = src_term.arrays();
    auto const& field_arrs = m_act_src(lev).const_arrays();
    amrex::ParallelFor(
        src_term, amrex::IntVect(0), AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
            src_arrs[nbx](i, j, k, n) += field_arrs[nbx](i, j, k, n);
        });
}

} // namespace kynema_sgf::pde::icns
