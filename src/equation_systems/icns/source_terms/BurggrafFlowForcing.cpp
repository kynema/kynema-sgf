#include "amr-wind/equation_systems/icns/source_terms/BurggrafFlowForcing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/physics/BurggrafFlow.H"

#include "AMReX_Gpu.H"

namespace amr_wind::pde::icns {

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
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto varr = m_bf_src(lev).const_array(mfi);

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        src_term(i, j, k, 0) += varr(i, j, k, 0);
        src_term(i, j, k, 1) += varr(i, j, k, 1);
        src_term(i, j, k, 2) += varr(i, j, k, 2);
    });
}

} // namespace amr_wind::pde::icns
