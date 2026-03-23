#ifndef FIELDUTILS_H
#define FIELDUTILS_H

#include "amr-wind/core/FieldDescTypes.H"
#include "AMReX_MultiFab.H"
#include "AMReX_Interpolater.H"

namespace amr_wind::field_impl {

/** Convert amr_wind::FieldLoc to index type for use with AMReX objects.
 *  \ingroup field_ops
 */
AMREX_INLINE amrex::IndexType index_type(const FieldLoc floc)
{
    switch (floc) {
    case FieldLoc::CELL:
        return amrex::IndexType::TheCellType();

    case FieldLoc::NODE:
        return amrex::IndexType::TheNodeType();

    case FieldLoc::XFACE:
        return amrex::IndexType(amrex::IntVect::TheDimensionVector(0));

    case FieldLoc::YFACE:
        return amrex::IndexType(amrex::IntVect::TheDimensionVector(1));

    case FieldLoc::ZFACE:
        return amrex::IndexType(amrex::IntVect::TheDimensionVector(2));
    }

    // Suppress warnings when compiling with CUDA
    return amrex::IndexType::TheCellType();
}

/** Given a base field name and a state, return a unique string that identifies
 *  this state in the field repository.
 *  \ingroup field_ops
 */
AMREX_INLINE std::string
field_name_with_state(const std::string& fname, const FieldState fstate)
{
    switch (fstate) {
    case FieldState::NP1:
        return fname;

    case FieldState::N:
        return fname + "__FS_Old";

    case FieldState::NM1:
        return fname + "__FS_NM1";

    case FieldState::NPH:
        return fname + "__FS_NPH";

    case FieldState::NMH:
        return fname + "__FS_NMH";
    }

    // Suppress warnings when compiling with CUDA
    return fname;
}

/** Check if the field is not a reserved field name
 *
 *  In amr_wind::FieldRepo, field with states are stored using suffixes that
 *  create a unique string to identify different states. This function checks if
 *  the user specified field name matches any of those reserved suffixes.
 *
 *. \return True if the field is not a reserved field name
 *  \ingroup field_ops
 */
AMREX_INLINE bool is_valid_field_name(const std::string& name)
{
    // Else make sure that our guard is not in the last 8 characters
    auto found = name.rfind("__FS_");
    return found == std::string::npos;
}

/** Return an AMReX coarse-to-fine field interpolation object based on
 *  amr_wind::FieldInterpolator
 *  \ingroup field_ops
 */
AMREX_INLINE amrex::Interpolater*
get_interpolation_operator(const FieldInterpolator itype)
{
    switch (itype) {
    case FieldInterpolator::PiecewiseConstant:
        return &amrex::pc_interp;

    case FieldInterpolator::CellConsLinear:
        return &amrex::cell_cons_interp;

    case FieldInterpolator::NodeBilinear:
        return &amrex::node_bilinear_interp;

    case FieldInterpolator::FaceDivFree:
        return &amrex::face_divfree_interp;

    case FieldInterpolator::FaceLinear:
        return &amrex::face_linear_interp;

    default:
        amrex::Abort("Unrecognized FieldInterpolator type");
    }

    return &amrex::cell_cons_interp;
}

AMREX_INLINE FieldState dof_state(const FieldState fstate)
{
    return (fstate == FieldState::New) ? FieldState::New : FieldState::Old;
}

AMREX_INLINE FieldState phi_state(const FieldState fstate)
{
    return (fstate == FieldState::Old) ? FieldState::Old : FieldState::NPH;
}

} // namespace amr_wind::field_impl

#endif /* FIELDUTILS_H */
