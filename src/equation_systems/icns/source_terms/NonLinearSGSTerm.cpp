#include "src/equation_systems/icns/source_terms/NonLinearSGSTerm.H"
#include "src/CFDSim.H"
#include "src/turbulence/LES/Kosovic.H"

#include "AMReX_Gpu.H"

namespace kynema_sgf::pde::icns {

NonLinearSGSTerm::NonLinearSGSTerm(const CFDSim& sim)
    : m_divNij(sim.repo().get_field("divNij"))
{}

NonLinearSGSTerm::~NonLinearSGSTerm() = default;

void NonLinearSGSTerm::operator()(
    const int lev, const FieldState /*fstate*/, amrex::MultiFab& src_term) const
{
    auto const& src_arrs = src_term.arrays();
    auto const& field_arrs = m_divNij(lev).const_arrays();
    amrex::ParallelFor(
        src_term, amrex::IntVect(0), AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
            src_arrs[nbx](i, j, k, n) += field_arrs[nbx](i, j, k, n);
        });
}

} // namespace kynema_sgf::pde::icns
