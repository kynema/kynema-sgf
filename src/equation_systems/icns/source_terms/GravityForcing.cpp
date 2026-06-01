#include "src/equation_systems/icns/source_terms/GravityForcing.H"
#include "src/CFDSim.H"
#include "src/core/FieldUtils.H"

#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::icns {

/** Gravity Forcing source term
 *
 *  Reads in the following parameters from `incflo` namespace:
 *
 *  - `gravity` acceleration due to gravity (m/s)
 */
GravityForcing::GravityForcing(const CFDSim& sim)
{
    amrex::ParmParse pp("incflo");
    pp.queryarr("gravity", m_gravity);
    pp.query("density", m_rho0_const);

    // Get density fields
    m_rho = &(sim.repo().get_field("density"));

    // Check if perturbational pressure desired
    amrex::ParmParse pp_icns("ICNS");
    pp_icns.query("use_perturb_pressure", m_use_perturb_pressure);
    m_use_reference_density = sim.repo().field_exists("reference_density");
    m_rho0 = m_use_reference_density
                 ? &(sim.repo().get_field("reference_density"))
                 : nullptr;
}

GravityForcing::~GravityForcing() = default;

/** Add the Gravity source term to the forcing array
 *
 *  \param lev AMR level
 *  \param fstate FieldState field
 *  \param src_term Cumulative forcing array
 */
void GravityForcing::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};

    auto const& src_arrs = src_term.arrays();
    auto const& rho_arrs =
        ((*m_rho).state(field_impl::phi_state(fstate)))(lev).const_arrays();
    auto const& rho0_arrs = m_use_reference_density
                                ? (*m_rho0)(lev).const_arrays()
                                : amrex::MultiArray4<amrex::Real const>();
    const bool ir0 = m_use_reference_density;
    const bool ipt = m_use_perturb_pressure;
    const amrex::Real mr0c = m_rho0_const;

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
            const amrex::Real factor =
                (!ipt ? 1.0_rt
                      : 1.0_rt - ((ir0 ? rho0_arrs[nbx](i, j, k) : mr0c) /
                                  rho_arrs[nbx](i, j, k)));
            src_arrs[nbx](i, j, k, n) += gravity[n] * factor;
        });
}

} // namespace kynema_sgf::pde::icns
