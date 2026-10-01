#ifndef HYBRID_LENGTH_SCALE_H
#define HYBRID_LENGTH_SCALE_H

#include <cmath>
#include "AMReX_GpuQualifiers.H"
#include "AMReX_Extension.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

/** Near-wall hybrid RANS-LES length scale
 *
 *  The grid-based LES length scale shrinks with every refinement level, while
 *  a wall-distance (RANS) length scale does not depend on the grid. Blending
 *  the two near the wall gives a length scale that is nearly independent of
 *  the refinement level close to the surface and reverts to the LES length
 *  scale aloft. These functions hold the blending used by the Kosovic model
 *  (surfaceRANS option) so that it can be shared with other LES models.
 */
namespace kynema_sgf::turbulence::hybrid_length {

/** Weight of the RANS length scale
 *
 *  \param height Height above the wall
 *  \param switch_height Height scale over which the weight decays
 *  \return exp(-height / switch_height), 1 at the wall and 0 far above it
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
rans_weight(const amrex::Real height, const amrex::Real switch_height)
{
    return std::exp(-height / switch_height);
}

/** Blended square of the length scale
 *
 *  \param les_length_sqr Square of the LES (grid-based) length scale
 *  \param rans_length_sqr Square of the RANS (wall-distance) length scale
 *  \param weight Weight of the RANS length scale from rans_weight
 *  \param exponent Exponent applied to the weights
 *  \return (1 - weight)^exponent * les_length_sqr + weight^exponent *
 *  rans_length_sqr
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real blended_length_sqr(
    const amrex::Real les_length_sqr,
    const amrex::Real rans_length_sqr,
    const amrex::Real weight,
    const amrex::Real exponent)
{
    return (std::pow(1.0_rt - weight, exponent) * les_length_sqr) +
           (std::pow(weight, exponent) * rans_length_sqr);
}

/** RANS length scale for a one-equation (k) subgrid model
 *
 *  For an eddy viscosity Ce * l * sqrt(k) with dissipation Ceps * k^(3/2) / l,
 *  production equals dissipation in the log layer when the eddy viscosity is
 *  kappa * u_* * z, which requires
 *  l = (kappa * z / phi_m) * (Ceps / Ce^3)^(1/4).
 *
 *  \param kappa von Karman constant
 *  \param height Height above the wall
 *  \param phi_m Monin-Obukhov stability function for momentum (1 if neutral)
 *  \param Ce Eddy viscosity coefficient of the one-equation model
 *  \param Ceps Dissipation coefficient of the one-equation model
 *  \return Length scale that recovers the log law in equilibrium
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real one_eq_rans_length(
    const amrex::Real kappa,
    const amrex::Real height,
    const amrex::Real phi_m,
    const amrex::Real Ce,
    const amrex::Real Ceps)
{
    return (kappa * height / phi_m) *
           std::sqrt(std::sqrt(Ceps / (Ce * Ce * Ce)));
}

} // namespace kynema_sgf::turbulence::hybrid_length

#endif /* HYBRID_LENGTH_SCALE_H */
