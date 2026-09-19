#ifndef IMMERSEDWALLMODEL_H
#define IMMERSEDWALLMODEL_H

#include "AMReX_GpuQualifiers.H"
#include "AMReX_REAL.H"

#include <cmath>

using namespace amrex::literals;

/** Helpers shared by the ImmersedTerrain physics and the immersed terrain
 *  source terms: the drag weight of a cell, the neighbor fraction that makes a
 *  face a wall, and the exact integration of a linear relaxation.
 */
namespace kynema_sgf::immersed_wall {

//! Weight of the immersed drag in a cell of terrain fraction beta.
//! ``fraction``: beta itself. ``center``: 1 when the cell center is inside the
//! terrain (beta >= threshold), else 0, so that mostly fluid partial cells are
//! treated as fluid cells and are not damped by the drag.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real solid_weight(
    const amrex::Real beta,
    const bool center_based,
    const amrex::Real threshold)
{
    if (center_based) {
        return (beta >= threshold) ? 1.0_rt : 0.0_rt;
    }
    return beta;
}

//! Exact integration of du/dt = -C (u - u_t) over one step: the source
//! -C_eff (u^n - u_t) with C_eff = (1 - exp(-C dt)) / dt reproduces
//! u^{n+1} = u_t + (u^n - u_t) exp(-C dt)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
exact_relaxation_rate(const amrex::Real C, const amrex::Real dt)
{
    return (1.0_rt - std::exp(-C * dt)) / dt;
}

/** Neighbor terrain fraction from which a face is a wall.
 *
 *  With center weighting a cell that is at least ``solid_threshold`` terrain
 *  is treated as solid, so a face next to it is a wall. With fraction
 *  weighting the partial cell carries its own wall, so only a neighbor that
 *  is entirely terrain makes a face a wall; a lower threshold would model the
 *  same wall again in the cell above the partial cell.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
wall_threshold(const bool center_based, const amrex::Real solid_threshold)
{
    return center_based ? solid_threshold : 1.0_rt;
}

} // namespace kynema_sgf::immersed_wall

#endif
