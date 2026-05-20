//
//  DensityBuoyancy.cpp
//  kynema-sgf
//

#include "src/equation_systems/icns/source_terms/DensityBuoyancy.H"
#include "src/CFDSim.H"
#include "src/core/FieldUtils.H"

#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::icns {

/** Density based buoyancy source term
 *
 *  Reads in the following parameters from `incflo` namespace:
 *
 *  - `gravity` acceleration due to gravity (m/s)
 *  - `reference density` Optional, default = `1.0_rt`
 */
DensityBuoyancy::DensityBuoyancy(const CFDSim& sim)
    : m_density(sim.repo().get_field("density"))
{
    // gravity in `incflo` namespace
    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
        pp.query("density", m_rho_0);
    }
}

DensityBuoyancy::~DensityBuoyancy() = default;

/** Add the Boussinesq source term to the forcing array
 *
 *  \param lev AMR level
 *  \param fstate field state
 *  \param src_term Cumulative forcing array
 */
void DensityBuoyancy::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const amrex::Real density_0 = m_rho_0;
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};

    FieldState den_state = field_impl::phi_state(fstate);
    auto const& src_arrs = src_term.arrays();
    auto const& density_arrs = m_density.state(den_state)(lev).const_arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
            const amrex::Real fac =
                1.0_rt - (density_0 / density_arrs[nbx](i, j, k));
            src_arrs[nbx](i, j, k, n) += gravity[n] * fac;
        });
}

} // namespace kynema_sgf::pde::icns
