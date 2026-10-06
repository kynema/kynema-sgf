#ifndef FORESTCANOPY_H
#define FORESTCANOPY_H

#include "AMReX_MultiFab.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::forestdrag {

/** Canopy turbulence inputs of the ForestDrag section
 *
 *  Read by the KLAxell turbulence model and its KransAxell TKE source. The
 *  canopy terms act where the forest_drag field (C_d * LAD) is positive, so
 *  they only need ForestDrag in canopy mode.
 */
struct CanopyTurbulence
{
    //! ForestDrag.canopy_tke: add the canopy TKE terms (off by default)
    bool m_enabled{false};
    //! ForestDrag.canopy_beta_p: wake production coefficient
    amrex::Real m_beta_p{1.0_rt};
    //! ForestDrag.canopy_beta_d: short-circuit dissipation coefficient
    amrex::Real m_beta_d{4.0_rt};
    //! ForestDrag.canopy_length_alpha: KLAxell length scale limit
    //! l <= alpha L_c; Harman & Finnigan (2007) give alpha = 2 beta^3 with
    //! beta = u*/U_h; 0.04 (beta = 0.27) best fits the Shaw & Schumann (1992)
    //! LES for LAI 2 and 5; zero disables the limit
    amrex::Real m_length_alpha{0.04_rt};
};

//! Read the canopy turbulence inputs from the ForestDrag section
CanopyTurbulence parse_canopy_turbulence();

/** Add the canopy TKE source to a TKE source term
 *
 *  S_k = f (beta_p |U|^3 - beta_d |U| k) with f = C_d LAD. The sink is
 *  integrated exactly over the time step, -k (1 - exp(-beta_d f |U| dt)) / dt,
 *  so it removes at most the local k in one step. Only cells with f > 0 are
 *  changed.
 *
 *  \param src_term TKE source term, valid cells updated
 *  \param forest_drag Forest drag field f = C_d LAD (1/m)
 *  \param velocity Velocity field
 *  \param tke Turbulent kinetic energy field
 *  \param beta_p Wake production coefficient
 *  \param beta_d Short-circuit dissipation coefficient
 *  \param dt Time step
 */
void add_canopy_tke_source(
    amrex::MultiFab& src_term,
    const amrex::MultiFab& forest_drag,
    const amrex::MultiFab& velocity,
    const amrex::MultiFab& tke,
    amrex::Real beta_p,
    amrex::Real beta_d,
    amrex::Real dt);

/** Limit the turbulence length scale inside the canopy
 *
 *  l <= alpha L_c with the canopy drag length L_c = 1 / f. The eddy viscosity,
 *  shear production and buoyancy production of KLAxell are linear in l, so
 *  they are rescaled by the same factor. Only cells with f > 0 and
 *  alpha L_c < l are changed.
 *
 *  \param forest_drag Forest drag field f = C_d LAD (1/m)
 *  \param alpha Length limit coefficient
 *  \param turb_lscale Turbulence length scale, valid cells updated
 *  \param mu_turb Turbulent viscosity, valid cells updated
 *  \param shear_prod Shear production, valid cells updated
 *  \param buoy_prod Buoyancy production, valid cells updated
 */
void limit_canopy_length_scale(
    const amrex::MultiFab& forest_drag,
    amrex::Real alpha,
    amrex::MultiFab& turb_lscale,
    amrex::MultiFab& mu_turb,
    amrex::MultiFab& shear_prod,
    amrex::MultiFab& buoy_prod);

} // namespace kynema_sgf::forestdrag

#endif
