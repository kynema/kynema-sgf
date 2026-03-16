#ifndef TRIG_OPS_H
#define TRIG_OPS_H

/** \file trig_ops.H
 *
 *  Trigonometric utilities
 */

#include <cmath>
#include <numbers>
#include "AMReX_REAL.H"
#include "AMReX_Gpu.H"

using namespace amrex::literals;

namespace amr_wind::utils {

//! Return \f$2 \pi\f$
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr amrex::Real two_pi()
{
    return 2.0_rt * std::numbers::pi_v<amrex::Real>;
}

//! Return \f$\pi / 2\f$
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr amrex::Real half_pi()
{
    return 0.5_rt * std::numbers::pi_v<amrex::Real>;
}

//! Convert from degrees to radians
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
radians(const amrex::Real deg_val)
{
    return std::numbers::pi_v<amrex::Real> * deg_val / 180.0_rt;
}

//! Convert from radians to degrees
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
degrees(const amrex::Real rad_val)
{
    return 180.0_rt * rad_val / std::numbers::pi_v<amrex::Real>;
}

} // namespace amr_wind::utils

#endif /* TRIG_OPS_H */
