#ifndef PASSIVESCALAR_H
#define PASSIVESCALAR_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/PDE.H"
#include "amr-wind/equation_systems/SourceTerm.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

/**
 *  \defgroup passive_eqn Passive scalar transport equation
 *  %Passive scalar transport equation
 *
 *  \ingroup eqsys
 */

/**
 *  \defgroup passive_src Passive Scalar Source terms
 *
 *  Source terms used with passive scalar transport equation
 *
 *  \ingroup passive_eqn
 */

/** Characteristics of passive scalar transport equation
 *  \ingroup passive_eqn
 */
struct PassiveScalar : ScalarTransport
{
    using MLDiffOp = amrex::MLABecLaplacian;
    using SrcTerm = SourceTerm;

    static std::string pde_name() { return "PassiveScalar"; }
    static std::string var_name() { return "passive_scalar"; }

    static constexpr amrex::Real default_bc_value = 0.0_rt;

    static constexpr int ndim = 1;
    static constexpr bool multiply_rho = true;
    static constexpr bool has_diffusion = true;
    static constexpr bool need_nph_state = true;
};

} // namespace amr_wind::pde

#endif
