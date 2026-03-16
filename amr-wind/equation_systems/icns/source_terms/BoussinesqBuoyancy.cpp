#include "amr-wind/equation_systems/icns/source_terms/BoussinesqBuoyancy.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldUtils.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::pde::icns {

/** Boussinesq buoyancy source term for ABL simulations
 *
 */
BoussinesqBuoyancy::BoussinesqBuoyancy(const CFDSim& sim)
    : m_temperature(sim.repo().get_field("temperature"))
    , m_transport(sim.transport_model())
{
    // gravity in `incflo` namespace
    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.queryarr("gravity", m_gravity);
}

BoussinesqBuoyancy::~BoussinesqBuoyancy() = default;

void BoussinesqBuoyancy::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};

    const auto& temp =
        m_temperature.state(field_impl::phi_state(fstate))(lev).const_array(
            mfi);
    amrex::FArrayBox beta_fab(bx, 1, amrex::The_Async_Arena());
    amrex::Array4<amrex::Real> const& beta_arr = beta_fab.array();
    m_transport.beta_impl(lev, mfi, bx, beta_arr);
    amrex::FArrayBox ref_theta_fab(bx, 1, amrex::The_Async_Arena());
    amrex::Array4<amrex::Real> const& ref_theta_arr = ref_theta_fab.array();
    m_transport.ref_theta_impl(lev, mfi, bx, ref_theta_arr);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real T = temp(i, j, k, 0);
        const amrex::Real T0 = ref_theta_arr(i, j, k);
        const amrex::Real fac = beta_arr(i, j, k) * (T0 - T);

        src_term(i, j, k, 0) += gravity[0] * fac;
        src_term(i, j, k, 1) += gravity[1] * fac;
        src_term(i, j, k, 2) += gravity[2] * fac;
    });
}

} // namespace amr_wind::pde::icns
