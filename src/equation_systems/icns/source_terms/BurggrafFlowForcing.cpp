#include "src/equation_systems/icns/source_terms/BurggrafFlowForcing.H"
#include "src/CFDSim.H"
#include "src/physics/BurggrafFlow.H"

#include "AMReX_Gpu.H"

namespace kynema_sgf::pde::icns {

BurggrafFlowForcing::BurggrafFlowForcing(const CFDSim& sim)
    : m_bf_src(sim.repo().get_field("bf_src_term"))
{
    if (!sim.physics_manager().contains("BurggrafFlow")) {
        amrex::Abort(
            "BurggrafFlowForcing requires BurggrafFlow physics to be active");
    }
}

BurggrafFlowForcing::~BurggrafFlowForcing() = default;

void BurggrafFlowForcing::operator()(
    const int lev, const FieldState /*fstate*/, amrex::MultiFab& src_term) const
{
    auto const& src_arrs = src_term.arrays();
    auto const& field_arrs = m_bf_src(lev).const_arrays();
    amrex::ParallelFor(
        src_term, amrex::IntVect(0), AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
            src_arrs[nbx](i, j, k, n) += field_arrs[nbx](i, j, k, n);
        });
}

} // namespace kynema_sgf::pde::icns
