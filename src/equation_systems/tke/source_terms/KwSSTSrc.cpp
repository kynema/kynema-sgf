#include <AMReX_Orientation.H>

#include "src/equation_systems/tke/source_terms/KwSSTSrc.H"
#include "src/CFDSim.H"
#include "src/turbulence/TurbulenceModel.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::tke {

KwSSTSrc::KwSSTSrc(const CFDSim& sim)
    : m_shear_prod(sim.repo().get_field("shear_prod"))
    , m_diss(sim.repo().get_field("dissipation"))
    , m_buoy_term(sim.repo().get_field("buoyancy_term"))
{}

KwSSTSrc::~KwSSTSrc() = default;

void KwSSTSrc::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const auto& shear_prod_arrs = (this->m_shear_prod)(lev).const_arrays();
    const auto& diss_arrs = (this->m_diss)(lev).const_arrays();
    const auto& buoy_arrs = (this->m_buoy_term)(lev).const_arrays();

    const amrex::Real factor = (fstate == FieldState::NPH) ? 0.5_rt : 1.0_rt;

    auto const& src_arrs = src_term.arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), 1,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
            src_arrs[nbx](i, j, k) += shear_prod_arrs[nbx](i, j, k) +
                                      (factor * diss_arrs[nbx](i, j, k)) +
                                      buoy_arrs[nbx](i, j, k);
        });
}

} // namespace kynema_sgf::pde::tke
