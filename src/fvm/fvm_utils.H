#ifndef FVM_UTILS_H
#define FVM_UTILS_H

#include "src/core/FieldRepo.H"
#include "src/fvm/stencils.H"
#include "src/utilities/math_ops.H"

namespace kynema_sgf::fvm::impl {

/** Apply a finite volume operator for a given field.
 *
 *  For each level, dispatches one fused kernel per stencil type over all FABs
 *  in the MultiFab, replacing the previous per-FAB MFIter loop.  The interior
 *  stencil is applied to every valid cell; boundary stencils overwrite only
 *  the cells that sit on a non-periodic domain face / edge / corner, with
 *  stencil applicability determined by each Stencil type's `applies_to()`
 *  method inside the kernel.
 */
template <typename FvmOp, typename FType>
AMREX_INLINE void apply(const FvmOp& fvmop, const FType& fld)
{
    namespace stencil = kynema_sgf::fvm::stencil;
    const int nlevels = fld.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = fld.repo().mesh().Geom(lev);
        const bool px = geom.isPeriodic(0);
        const bool py = geom.isPeriodic(1);
        const bool pz = geom.isPeriodic(2);

        fvmop.template applyFused<stencil::StencilInterior>(lev);

        // faces
        if (!px) {
            fvmop.template applyFused<stencil::StencilILO>(lev);
            fvmop.template applyFused<stencil::StencilIHI>(lev);
        }
        if (!py) {
            fvmop.template applyFused<stencil::StencilJLO>(lev);
            fvmop.template applyFused<stencil::StencilJHI>(lev);
        }
        if (!pz) {
            fvmop.template applyFused<stencil::StencilKLO>(lev);
            fvmop.template applyFused<stencil::StencilKHI>(lev);
        }

        // edges
        if (!px && !py) {
            fvmop.template applyFused<stencil::StencilIhiJlo>(lev);
            fvmop.template applyFused<stencil::StencilIhiJhi>(lev);
            fvmop.template applyFused<stencil::StencilIloJlo>(lev);
            fvmop.template applyFused<stencil::StencilIloJhi>(lev);
        }
        if (!px && !pz) {
            fvmop.template applyFused<stencil::StencilIhiKlo>(lev);
            fvmop.template applyFused<stencil::StencilIhiKhi>(lev);
            fvmop.template applyFused<stencil::StencilIloKlo>(lev);
            fvmop.template applyFused<stencil::StencilIloKhi>(lev);
        }
        if (!py && !pz) {
            fvmop.template applyFused<stencil::StencilJhiKlo>(lev);
            fvmop.template applyFused<stencil::StencilJhiKhi>(lev);
            fvmop.template applyFused<stencil::StencilJloKlo>(lev);
            fvmop.template applyFused<stencil::StencilJloKhi>(lev);
        }

        // corners
        if (!px && !py && !pz) {
            fvmop.template applyFused<stencil::StencilIloJloKlo>(lev);
            fvmop.template applyFused<stencil::StencilIloJloKhi>(lev);
            fvmop.template applyFused<stencil::StencilIloJhiKlo>(lev);
            fvmop.template applyFused<stencil::StencilIloJhiKhi>(lev);
            fvmop.template applyFused<stencil::StencilIhiJloKlo>(lev);
            fvmop.template applyFused<stencil::StencilIhiJloKhi>(lev);
            fvmop.template applyFused<stencil::StencilIhiJhiKlo>(lev);
            fvmop.template applyFused<stencil::StencilIhiJhiKhi>(lev);
        }
    }
}

} // namespace kynema_sgf::fvm::impl

#endif /* FVM_UTILS_H */
