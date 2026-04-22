#ifndef TKE_H
#define TKE_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/PDE.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/equation_systems/tke/TKESource.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

/**
 *  \defgroup tke_eqn Turbulent Kinetic Energy (TKE) transport equation
 *  %TKE transport equation
 *
 *  \ingroup eqsys
 */

/**
 *  \defgroup tke_src TKE Source terms
 *
 *  Source terms used with TKE transport equation
 *
 *  \ingroup tke_eqn
 */

/** Characteristics of TKE transport equation
 *
 *  \ingroup tke_eqn
 */
struct TKE : ScalarTransport
{
    using MLDiffOp = amrex::MLABecLaplacian;
    using SrcTerm = TKESource;

    static std::string pde_name() { return "TKE"; }
    static std::string var_name() { return "tke"; }

    static constexpr amrex::Real default_bc_value = 0.0_rt;

    static constexpr int ndim = 1;
    static constexpr bool multiply_rho = true;
    static constexpr bool has_diffusion = true;
    static constexpr bool need_nph_state = true;
};

} // namespace amr_wind::pde

#endif /* TKE_H */
