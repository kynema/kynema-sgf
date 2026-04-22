#ifndef WAVES_UTILS_H
#define WAVES_UTILS_H

#include <cmath>
#include <numbers>
#include "amr-wind/utilities/constants.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::ocean_waves::utils {

using WaveVec = amrex::GpuArray<amrex::Real, 4>;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real free_surface_to_vof(
    const amrex::Real eta, const amrex::Real z, const amrex::Real dz)
{
    amrex::Real volfrac = 0.0_rt;
    if (eta - z >= dz / 2.0_rt) {
        volfrac = 1.0_rt;
    } else if (eta - z <= -dz / 2.0_rt) {
        volfrac = 0.0_rt;
    } else if (std::abs(eta - z) < dz / 2.0_rt) {
        if (eta <= z) {
            volfrac = 1.0_rt - ((z - eta + dz / 2.0_rt) / dz);
        } else {
            volfrac = (eta - z + dz / 2.0_rt) / dz;
        }
    }

    return volfrac;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
gamma_generate(const amrex::Real x, const amrex::Real gen_length)
{
    const amrex::Real xtilde = amrex::max<amrex::Real>(
        amrex::min<amrex::Real>(1.0_rt - (x / gen_length), 1.0_rt), 0.0_rt);
    return (
        1.0_rt - (std::expm1(std::pow(xtilde, 3.5_rt)) / std::expm1(1.0_rt)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real gamma_absorb(
    const amrex::Real x,
    const amrex::Real absorb_length,
    const amrex::Real absorb_length_factor)
{
    const amrex::Real xtilde = amrex::max<amrex::Real>(
        amrex::min<amrex::Real>(
            x / (absorb_length * absorb_length_factor), 1.0_rt),
        0.0_rt);
    return (
        1.0_rt - (std::expm1(std::pow(xtilde, 3.5_rt)) / std::expm1(1.0_rt)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
ramp(const amrex::Real time, const amrex::Real ramp_period)
{
    amrex::Real f = 1.0_rt;
    if (time < ramp_period) {
        f = (time / ramp_period) -
            ((1.0_rt / std::numbers::pi_v<amrex::Real>)*std::sin(
                std::numbers::pi_v<amrex::Real> * (time / ramp_period)));
    }
    return f;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real combine_linear(
    const amrex::Real factor,
    const amrex::Real target,
    const amrex::Real current)
{
    return ((1.0_rt - factor) * target) + (factor * current);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE WaveVec harmonize_profiles_1d(
    const amrex::Real x,
    const amrex::Real left_bdy,
    const amrex::Real left_length,
    const amrex::Real right_bdy,
    const amrex::Real right_length,
    const WaveVec left,
    const WaveVec bulk,
    const WaveVec right)
{
    // Initialize boundary regions entirely with boundary solutions
    if (x <= left_bdy + left_length) {
        return left;
    }
    if (x + right_length >= right_bdy) {
        return right;
    }
    // Combine solutions in central region using relax zone length scales
    const amrex::Real Gamma_left =
        gamma_generate(x - (left_bdy + left_length), 0.5_rt * left_length);
    const amrex::Real Gamma_right = gamma_absorb(
        x - (right_bdy - 1.5_rt * right_length), 0.5_rt * right_length, 1.0_rt);
    const bool right_inactive = Gamma_right > 1.0_rt - constants::TIGHT_TOL;
    WaveVec combo;
    for (int n = 0; n < 4; ++n) {
        combo[n] = right_inactive
                       ? combine_linear(Gamma_left, left[n], bulk[n])
                       : combine_linear(Gamma_right, right[n], bulk[n]);
    }
    return combo;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE WaveVec harmonize_profiles_2d(
    const amrex::Real x,
    const amrex::Real y,
    const amrex::Real left_bdy,
    const amrex::Real left_length,
    const amrex::Real right_bdy,
    const amrex::Real right_length,
    const amrex::Real lo_y_bdy,
    const amrex::Real hi_y_bdy,
    const amrex::Real y_length,
    const WaveVec left,
    const WaveVec bulk,
    const WaveVec right)
{
    // First get forcing profile, applies in x and y zones
    auto profile_forcing = harmonize_profiles_1d(
        x, left_bdy, left_length, right_bdy, right_length, left, left, right);

    // Also get initial profile according to x direction
    auto profile_init = harmonize_profiles_1d(
        x, left_bdy, left_length, right_bdy, right_length, left, bulk, right);

    // Then harmonize in y direction if applicable
    if (y_length >= constants::EPS) {
        profile_init = harmonize_profiles_1d(
            y, lo_y_bdy, y_length, hi_y_bdy, y_length, profile_forcing,
            profile_init, profile_forcing);
    }

    return profile_init;
}

} // namespace amr_wind::ocean_waves::utils

#endif
