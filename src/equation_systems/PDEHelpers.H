#ifndef PDEHELPERS_H
#define PDEHELPERS_H

#include <string>

#include "amr-wind/equation_systems/PDEFields.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/FieldFillPatchOps.H"
#include "amr-wind/core/FieldBCOps.H"

namespace amr_wind {
namespace pde_impl {

//! Field name that holds the convective term
AMREX_INLINE std::string conv_term_name(const std::string& var)
{
    return var + "_conv_term";
}

//! Field name that holds the diffusion term (divtau and laplacian for scalars)
AMREX_INLINE std::string diff_term_name(const std::string& var)
{
    return var + "_diff_term";
}

//! Field name for all the source terms
AMREX_INLINE std::string src_term_name(const std::string& var)
{
    return var + "_src_term";
}

//! Effective viscosity for the transport equation
AMREX_INLINE std::string mueff_name(const std::string& var)
{
    return var + "_mueff";
}

} // namespace pde_impl

namespace pde {

/** Declare fields necessary for solving a PDE
 *
 *  This function defines the boilerplate code to declare all the fields in
 *  PDEFields that are necessary to solve the equation system.
 *
 *  \tparam PDE Trait defining this equation system (e.g., ICNS, Temperature)
 *  \tparam Scheme The numerical scheme (e.g., Godnuov, MOL)
 */
template <typename PDE, typename Scheme>
PDEFields create_fields_instance(
    const SimTime& time,
    FieldRepo& repo,
    const FieldInterpolator itype = FieldInterpolator::CellConsLinear)
{
    repo.declare_field(
        PDE::var_name(), PDE::ndim, Scheme::nghost_state, Scheme::num_states);
    repo.declare_field(pde_impl::mueff_name(PDE::var_name()), 1, 1, 1);
    repo.declare_field(
        pde_impl::src_term_name(PDE::var_name()), PDE::ndim, Scheme::nghost_src,
        1);
    repo.declare_field(
        pde_impl::diff_term_name(PDE::var_name()), PDE::ndim, 0,
        Scheme::num_diff_states);
    repo.declare_field(
        pde_impl::conv_term_name(PDE::var_name()), PDE::ndim, 0,
        Scheme::num_conv_states);

    PDEFields fields(repo, PDE::var_name());
    fields.field.register_fill_patch_op<FieldFillPatchOps<FieldBCDirichlet>>(
        repo.mesh(), time, itype);
    fields.src_term.register_fill_patch_op<FieldFillPatchOps<FieldBCNoOp>>(
        repo.mesh(), time, FieldInterpolator::PiecewiseConstant);
    fields.mueff.set_default_fillpatch_bc(time);

    fields.field.fillpatch_on_regrid() = true;

    if (PDE::need_nph_state) {
        fields.field.create_state(FieldState::NPH);
    }

    return fields;
}

} // namespace pde
} // namespace amr_wind

#endif /* PDEHELPERS_H */
