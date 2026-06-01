#ifndef SDR_H
#define SDR_H

#include "src/equation_systems/PDETraits.H"
#include "src/equation_systems/SchemeTraits.H"
#include "src/equation_systems/PDEHelpers.H"
#include "src/equation_systems/PDE.H"
#include "src/core/field_ops.H"
#include "src/equation_systems/sdr/SDRSource.H"

namespace kynema_sgf::pde {

/**
 *  \defgroup sdr_eqn Specific Dissipation Rate (SDR) transport equation
 *  %SDR transport equation
 *
 *  \ingroup eqsys
 */

/**
 *  \defgroup sdr_src SDR source terms
 *
 *  Source terms used with SDR transport equation
 *
 *  \ingroup sdr_eqn
 */

/** Characteristics of SDR transport equation
 *  \ingroup sdr_eqn
 */
struct SDR : ScalarTransport
{
    using MLDiffOp = amrex::MLABecLaplacian;
    using SrcTerm = SDRSource;

    static std::string pde_name() { return "SDR"; }
    static std::string var_name() { return "sdr"; }

    static constexpr amrex::Real default_bc_value = 0.0_rt;

    static constexpr int ndim = 1;
    static constexpr bool multiply_rho = true;
    static constexpr bool has_diffusion = true;
    static constexpr bool need_nph_state = true;
};

} // namespace kynema_sgf::pde

#endif /* SDR_H */
