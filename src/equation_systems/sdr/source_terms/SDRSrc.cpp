#include <AMReX_Orientation.H>

#include "src/equation_systems/sdr/source_terms/SDRSrc.H"
#include "src/CFDSim.H"
#include "src/core/FieldUtils.H"
#include "src/turbulence/TurbulenceModel.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::tke {

SDRSrc::SDRSrc(const CFDSim& sim)
    : m_time(sim.time())
    , m_sdr_src(sim.repo().get_field("omega_src"))
    , m_sdr_diss(sim.repo().get_field("sdr_dissipation"))
    , m_sdr(sim.repo().get_field("sdr"))
    , m_lhs(sim.repo().get_field("sdr_lhs_src_term"))
    , m_density(sim.repo().get_field("density"))
{}

SDRSrc::~SDRSrc() = default;

void SDRSrc::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const auto& sdr_src_arrs = (this->m_sdr_src)(lev).const_arrays();
    const auto& sdr_diss_arrs = (this->m_sdr_diss)(lev).const_arrays();
    // The model terms are density weighted and the SDR source is multiplied by
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
                    (sdr_diss_arrs[nbx](i, j, k) + sdr_src_arrs[nbx](i, j, k)) /
                    rho_arrs[nbx](i, j, k);
            });
        return;
    }

    // Right-hand side of the SDR update. The turbulence model put
    // lhs_src_term * sdr^{n+1} on the left-hand side, so lhs_src_term * sdr is
    // added here at the state the model was evaluated on (delta form).
    // lhs_src_term is zero for explicit diffusion.
    const amrex::Real inv_dt = 1.0_rt / m_time.delta_t();
    const auto& sdr_arrs =
        m_sdr.state(field_impl::dof_state(fstate))(lev).const_arrays();
    const auto& lhs_arrs = m_lhs(lev).const_arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), 1,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
            src_arrs[nbx](i, j, k) +=
                (sdr_diss_arrs[nbx](i, j, k) + sdr_src_arrs[nbx](i, j, k) +
                 (lhs_arrs[nbx](i, j, k) * inv_dt * sdr_arrs[nbx](i, j, k))) /
                rho_arrs[nbx](i, j, k);
        });
}

} // namespace kynema_sgf::pde::tke
