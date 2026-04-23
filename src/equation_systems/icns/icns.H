#ifndef ICNS_H
#define ICNS_H

#include "src/equation_systems/PDETraits.H"
#include "src/equation_systems/SchemeTraits.H"
#include "src/equation_systems/PDEHelpers.H"
#include "src/equation_systems/PDE.H"
#include "src/equation_systems/icns/MomentumSource.H"

#include "AMReX_MLTensorOp.H"

namespace kynema_sgf::pde {

/**
 *  \defgroup icns Incompressible Navier-Stokes
 *
 *  Incompressible Navier-Stokes equation system
 *
 *  \ingroup eqsys
 */

/**
 *  \defgroup icns_src ICNS Source terms
 *
 *  Source terms used with incompressible Navier-Stokes momentum equation
 *
 *  \ingroup icns
 */

/** Characteristics of ICNS equation system
 *  \ingroup icns
 */
struct ICNS : VectorTransport
{
    using SrcTerm = MomentumSource;

    static std::string pde_name() { return "ICNS"; }
    static std::string var_name() { return "velocity"; }

    static constexpr int ndim = AMREX_SPACEDIM;

    static constexpr bool multiply_rho = true;
    static constexpr bool has_diffusion = true;

    static constexpr bool need_nph_state = true;
};

} // namespace kynema_sgf::pde

#endif /* ICNS_H */
