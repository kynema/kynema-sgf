#include "src/equation_systems/icns/source_terms/BoussinesqBuoyancy.H"
#include "src/CFDSim.H"
#include "src/core/FieldUtils.H"

#include "AMReX_ParmParse.H"

namespace kynema_sgf::pde::icns {

BoussinesqBuoyancy::BoussinesqBuoyancy(const CFDSim& sim)
    : m_temperature(sim.repo().get_field("temperature"))
    , m_repo(sim.repo())
    , m_transport(sim.transport_model())
{
    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.queryarr("gravity", m_gravity);
}

BoussinesqBuoyancy::~BoussinesqBuoyancy() = default;

void BoussinesqBuoyancy::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    if (!m_beta_scratch) {
        m_beta_scratch = m_repo.create_scratch_field(1, 0);
    }
    if (!m_ref_theta_scratch) {
        m_ref_theta_scratch = m_repo.create_scratch_field(1, 0);
    }

    m_transport.beta_fill(lev, (*m_beta_scratch)(lev));
    m_transport.ref_theta_fill(lev, (*m_ref_theta_scratch)(lev));

    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};

    auto const& src_arrs = src_term.arrays();
    auto const& temp_arrs =
        m_temperature.state(field_impl::phi_state(fstate))(lev).const_arrays();
    auto const& beta_arrs = (*m_beta_scratch)(lev).const_arrays();
    auto const& ref_arrs = (*m_ref_theta_scratch)(lev).const_arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
            const amrex::Real T = temp_arrs[nbx](i, j, k, 0);
            const amrex::Real T0 = ref_arrs[nbx](i, j, k);
            const amrex::Real fac = beta_arrs[nbx](i, j, k) * (T0 - T);
            src_arrs[nbx](i, j, k, n) += gravity[n] * fac;
        });
}

} // namespace kynema_sgf::pde::icns
