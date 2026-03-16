#ifndef TEMPERATURE_H
#define TEMPERATURE_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/PDE.H"
#include "amr-wind/equation_systems/temperature/TemperatureSource.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

/**
 *  \defgroup temp_eqn Temperature transport equation
 *  %Temperature transport equation
 *
 *  \ingroup eqsys
 */

/**
 *  \defgroup temp_src Temperature Source terms
 *
 *  Source terms used with temperature transport equation
 *
 *  \ingroup temp_eqn
 */

/** Characteristics of temperature transport equation
 *  \ingroup temp_eqn
 */
struct Temperature : ScalarTransport
{
    using MLDiffOp = amrex::MLABecLaplacian;
    using SrcTerm = TemperatureSource;

    static std::string pde_name() { return "Temperature"; }
    static std::string var_name() { return "temperature"; }

    static constexpr amrex::Real default_bc_value = 0.0_rt;

    static constexpr int ndim = 1;
    static constexpr bool multiply_rho = true;
    static constexpr bool has_diffusion = true;
    static constexpr bool need_nph_state = true;
};

/** Effective thermal diffusivity update operator
 *  \ingroup temp_eqn
 */
template <>
struct TurbulenceOp<Temperature>
{
    TurbulenceOp(turbulence::TurbulenceModel& tmodel, PDEFields& fields)
        : m_tmodel(tmodel), m_fields(fields)
    {}

    void operator()() { m_tmodel.update_alphaeff(m_fields.mueff); }

    turbulence::TurbulenceModel& m_tmodel;
    PDEFields& m_fields;
};

} // namespace amr_wind::pde

#endif /* TEMPERATURE_H */
