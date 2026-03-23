#include <AMReX_Orientation.H>

#include "amr-wind/equation_systems/tke/source_terms/KwSSTSrc.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::tke {

KwSSTSrc::KwSSTSrc(const CFDSim& sim)
    : m_shear_prod(sim.repo().get_field("shear_prod"))
    , m_diss(sim.repo().get_field("dissipation"))
    , m_buoy_term(sim.repo().get_field("buoyancy_term"))
{}

KwSSTSrc::~KwSSTSrc() = default;

void KwSSTSrc::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& shear_prod_arr = (this->m_shear_prod)(lev).array(mfi);
    const auto& diss_arr = (this->m_diss)(lev).array(mfi);
    const auto& buoy_arr = (this->m_buoy_term)(lev).array(mfi);

    const amrex::Real factor = (fstate == FieldState::NPH) ? 0.5_rt : 1.0_rt;
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        src_term(i, j, k) += shear_prod_arr(i, j, k) +
                             (factor * diss_arr(i, j, k)) + buoy_arr(i, j, k);
    });
}

} // namespace amr_wind::pde::tke
