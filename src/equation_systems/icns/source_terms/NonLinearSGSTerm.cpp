#include "amr-wind/equation_systems/icns/source_terms/NonLinearSGSTerm.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/turbulence/LES/Kosovic.H"

#include "AMReX_Gpu.H"

namespace amr_wind::pde::icns {

NonLinearSGSTerm::NonLinearSGSTerm(const CFDSim& sim)
    : m_divNij(sim.repo().get_field("divNij"))
{}

NonLinearSGSTerm::~NonLinearSGSTerm() = default;

void NonLinearSGSTerm::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto varr = m_divNij(lev).const_array(mfi);

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        src_term(i, j, k, 0) += varr(i, j, k, 0);
        src_term(i, j, k, 1) += varr(i, j, k, 1);
        src_term(i, j, k, 2) += varr(i, j, k, 2);
    });
}

} // namespace amr_wind::pde::icns
