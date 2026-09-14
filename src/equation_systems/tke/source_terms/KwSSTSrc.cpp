#include <AMReX_Orientation.H>

#include "src/equation_systems/tke/source_terms/KwSSTSrc.H"
#include "src/CFDSim.H"
#include "src/core/FieldUtils.H"
#include "src/turbulence/TurbulenceModel.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::tke {

KwSSTSrc::KwSSTSrc(const CFDSim& sim)
    : m_time(sim.time())
    , m_shear_prod(sim.repo().get_field("shear_prod"))
    , m_diss(sim.repo().get_field("dissipation"))
    , m_buoy_term(sim.repo().get_field("buoyancy_term"))
    , m_tke(sim.repo().get_field("tke"))
    , m_lhs(sim.repo().get_field("tke_lhs_src_term"))
    , m_density(sim.repo().get_field("density"))
{}

KwSSTSrc::~KwSSTSrc() = default;

void KwSSTSrc::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const auto& shear_prod_arrs = (this->m_shear_prod)(lev).const_arrays();
    const auto& diss_arrs = (this->m_diss)(lev).const_arrays();
    const auto& buoy_arrs = (this->m_buoy_term)(lev).const_arrays();
    // The model terms are density weighted and the TKE source is multiplied by
    // this density after this call, so the terms are divided by it here
    const auto& rho_arrs =
        m_density.state(field_impl::phi_state(fstate))(lev).const_arrays();

    auto const& src_arrs = src_term.arrays();

    if (fstate == FieldState::Old) {
        // Forcing for the Godunov predictor: the full explicit balance
        amrex::ParallelFor(
            src_term, amrex::IntVect(0), 1,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
                src_arrs[nbx](i, j, k) +=
                    (shear_prod_arrs[nbx](i, j, k) + diss_arrs[nbx](i, j, k) +
                     buoy_arrs[nbx](i, j, k)) /
                    rho_arrs[nbx](i, j, k);
            });
        return;
    }

    // Right-hand side of the TKE update. The turbulence model put
    // lhs_src_term * tke^{n+1} on the left-hand side, so lhs_src_term * tke is
    // added here at the state the model was evaluated on (delta form).
    // lhs_src_term is zero for explicit diffusion.
    const amrex::Real inv_dt = 1.0_rt / m_time.delta_t();
    const auto& tke_arrs =
        m_tke.state(field_impl::dof_state(fstate))(lev).const_arrays();
    const auto& lhs_arrs = m_lhs(lev).const_arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), 1,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
            src_arrs[nbx](i, j, k) +=
                (shear_prod_arrs[nbx](i, j, k) + diss_arrs[nbx](i, j, k) +
                 buoy_arrs[nbx](i, j, k) +
                 (lhs_arrs[nbx](i, j, k) * inv_dt * tke_arrs[nbx](i, j, k))) /
                rho_arrs[nbx](i, j, k);
        });
}

} // namespace kynema_sgf::pde::tke
