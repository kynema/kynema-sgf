#ifndef PDETRAITS_H
#define PDETRAITS_H

namespace amr_wind::pde {

/** Characteristics of a vector transport equation
 *  \ingroup eqsys
 */
struct VectorTransport
{
    // Linear operator for the diffusion solve
    // using MLDiffOp = amrex::MLABecLaplacian;

    // Number of components per degree-of-freedom
    // static constexpr int ndim = AMREX_SPACEDIM;

    // Flag indicating whether this equation will use `(n + 1/2)` state
    // static constexpr bool need_nph_state = true;

    // Flag indicating if the terms need to be multiplied by density (e.g.,
    // source term)
    // static constexpr bool multiply_rho = true;

    // Flag indicating whether the equation has a diffusion term
    // static constexpr bool has_diffusion = true;
};

/** Characteristics of a scalar transport equation
 *  \ingroup eqsys
 *
 *  This class defines common characteristics of a scalar transport equation
 *  that concrete traits can override. For example, compare with Density trait
 */
struct ScalarTransport
{
    // AMReX linear operator used to solve the diffusion equation
    // using MLDiffOp = amrex::MLABecLaplacian;
    // Base class of the source term used to create specific instances
    // using SrcTerm = SourceTerm;

    // Number of components for the PDE
    // static constexpr int ndim = 1;

    // Do the PDE terms have to be multiplied by density
    // static constexpr bool multiply_rho = true;

    // Does this PDE have a diffusion term
    // static constexpr bool has_diffusion = true;

    // Default value to initialize the BC array
    // static constexpr amrex::Real default_bc_value = 0.0_rt;

    // Does this scalar need an NPH state
    // static constexpr bool need_nph_state = true;
};

} // namespace amr_wind::pde

#endif /* PDETRAITS_H */
