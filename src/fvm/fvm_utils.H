#ifndef FVM_UTILS_H
#define FVM_UTILS_H

#include "src/core/FieldRepo.H"
#include "src/fvm/stencils.H"
#include "src/utilities/math_ops.H"

#include "AMReX_MFIter.H"

namespace kynema_sgf::fvm::impl {

/** Apply a finite volume operator for a given field.
 *
 *  Dispatches a single fused kernel per level over all FABs in the MultiFab,
 *  with stencil-coefficient selection (interior / face / edge / corner) handled
 *  inside the kernel based on each cell's position relative to the domain
 *  boundary.  This replaces the previous approach of up to 27 separate kernel
 *  launches per FAB.
 */
template <typename FvmOp, typename FType>
AMREX_INLINE void apply(const FvmOp& fvmop, const FType& fld)
{
    const int nlevels = fld.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        fvmop.applyFused(lev);
    }
}

} // namespace kynema_sgf::fvm::impl

#endif /* FVM_UTILS_H */
