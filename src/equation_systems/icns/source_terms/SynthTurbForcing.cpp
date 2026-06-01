#include "src/equation_systems/icns/source_terms/SynthTurbForcing.H"
#include "src/CFDSim.H"

#include "AMReX_Gpu.H"

namespace kynema_sgf::pde::icns {

/** Synthetic Turbulence forcing term
 *
 *
 *
 */
SynthTurbForcing::SynthTurbForcing(const CFDSim& sim)
    : m_turb_force(sim.repo().get_field("synth_turb_forcing"))
{
    if (!sim.physics_manager().contains("SyntheticTurbulence")) {
        amrex::Abort(
            "SynthTurbForcing: SyntheticTurbulence physics not enabled");
    }
}

SynthTurbForcing::~SynthTurbForcing() = default;

void SynthTurbForcing::operator()(
    const int lev, const FieldState /*fstate*/, amrex::MultiFab& src_term) const
{
    auto const& src_arrs = src_term.arrays();
    auto const& field_arrs = m_turb_force(lev).arrays();
    amrex::ParallelFor(
        src_term, amrex::IntVect(0), AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
            src_arrs[nbx](i, j, k, n) += field_arrs[nbx](i, j, k, n);
        });
}

} // namespace kynema_sgf::pde::icns
