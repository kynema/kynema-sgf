#ifndef LEVELSET_H
#define LEVELSET_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/PDE.H"
#include "amr-wind/equation_systems/SourceTerm.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

/**
 *  \defgroup levelset Levelset transport equation
 *  \ingroup eqsys
 */

/**
 *  \defgroup levelset_src Levelset source terms
 *
 *  Source terms used in Levelset transport equation
 *  \ingroup levelset
 */

/** Characteristics of Levelset transport equation system
 *  \ingroup levelset
 */
struct Levelset : ScalarTransport
{
    using MLDiffOp = amrex::MLABecLaplacian;
    using SrcTerm = SourceTerm;

    static std::string pde_name() { return "Levelset"; }
    static std::string var_name() { return "levelset"; }

    static constexpr int ndim = 1;

    static constexpr bool multiply_rho = false;
    static constexpr bool has_diffusion = false;
    static constexpr bool need_nph_state = true;

    static constexpr amrex::Real default_bc_value = 0.0_rt;
};

} // namespace amr_wind::pde

#endif /* LEVELSET_H */
