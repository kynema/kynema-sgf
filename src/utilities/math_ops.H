#ifndef MATH_OPS_H
#define MATH_OPS_H

/** \file math_ops.H
 *
 *  Basic math utilities
 */

#include "AMReX_REAL.H"
#include "AMReX_Gpu.H"

using namespace amrex::literals;

namespace kynema_sgf::utils {

//! Helper function to do pow() with integer exponents and output amrex::Real
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
powi(const amrex::Real& x, const int& y)
{
    return static_cast<amrex::Real>(std::pow(x, y));
}

//! Prevent other argument types from accidentally being used
template <typename T, typename U>
amrex::Real powi(T, U) = delete;

} // namespace kynema_sgf::utils

#endif /* MATH_OPS_H */
