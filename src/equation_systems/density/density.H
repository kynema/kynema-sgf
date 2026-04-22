#ifndef DENSITY_H
#define DENSITY_H

#include "amr-wind/equation_systems/PDETraits.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/equation_systems/PDE.H"
#include "amr-wind/equation_systems/density/DensitySource.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

struct Density : ScalarTransport
{
    using MLDiffOp = amrex::MLABecLaplacian;
    using SrcTerm = DensitySource;

    static std::string pde_name() { return "Density"; }
    static std::string var_name() { return "density"; }

    static constexpr int ndim = 1;

    static constexpr bool multiply_rho = false;
    static constexpr bool has_diffusion = false;
    static constexpr bool need_nph_state = true;

    static constexpr amrex::Real default_bc_value = 1.0_rt;
};

} // namespace amr_wind::pde

#endif /* DENSITY_H */
