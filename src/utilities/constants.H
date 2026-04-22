#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "AMReX_REAL.H"
#include "AMReX_Math.H"
#include "AMReX_GpuQualifiers.H"

using namespace amrex::literals;

namespace amr_wind::constants {

//! A number very close to zero
static constexpr amrex::Real EPS =
    static_cast<amrex::Real>(std::numeric_limits<amrex::Real>::epsilon());

//! A tight tolerance
static constexpr amrex::Real TIGHT_TOL =
    std::numeric_limits<amrex::Real>::epsilon() * 1.0e3_rt;

//! A loose tolerance
static constexpr amrex::Real LOOSE_TOL =
    std::numeric_limits<amrex::Real>::epsilon() * 1.0e8_rt;

//! A large negative number
static constexpr amrex::Real LOW_NUM = static_cast<amrex::Real>(
    std::numeric_limits<amrex::Real>::lowest() * 1.0e-10_rt);

//! A large positive number
static constexpr amrex::Real LARGE_NUM = static_cast<amrex::Real>(
    std::numeric_limits<amrex::Real>::max() * 1.0e-10_rt);

//! Closeness comparison
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
is_close(const amrex::Real a, const amrex::Real b)
{
    return amrex::Math::abs(a - b) < EPS;
}

//! Boltzmann constant (J/K)
static constexpr amrex::Real BOLTZMANN_CONSTANT = 1.380649_rt * 1.0e-23_rt;

//! Avogadro's constant (mol^-1)
static constexpr amrex::Real AVOGADRO_CONSTANT = 6.02214076_rt * 1.0e23_rt;

//! Universal gas constant (J/K mol)
static constexpr amrex::Real UNIVERSAL_GAS_CONSTANT =
    AVOGADRO_CONSTANT * BOLTZMANN_CONSTANT;

//! Heat capacity ratio
static constexpr amrex::Real HEAT_CAPACITY_RATIO = 1.4_rt;

//! Molar mass of air (kg/mol)
static constexpr amrex::Real MOLAR_MASS_AIR = 0.02896492_rt;

//! Molar mass of water vapor (kg/mol)
static constexpr amrex::Real MOLAR_MASS_WATER_VAPOR = 0.01801528_rt;

} // namespace amr_wind::constants
#endif
