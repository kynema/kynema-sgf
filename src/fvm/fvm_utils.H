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

        fvmop.template apply<stencil::StencilInterior>(lev);

        // faces
        if (!px) {
            fvmop.template apply<stencil::StencilILO>(lev);
            fvmop.template apply<stencil::StencilIHI>(lev);
        }
        if (!py) {
            fvmop.template apply<stencil::StencilJLO>(lev);
            fvmop.template apply<stencil::StencilJHI>(lev);
        }
        if (!pz) {
            fvmop.template apply<stencil::StencilKLO>(lev);
            fvmop.template apply<stencil::StencilKHI>(lev);
        }

        // edges
        if (!px && !py) {
            fvmop.template apply<stencil::StencilIhiJlo>(lev);
            fvmop.template apply<stencil::StencilIhiJhi>(lev);
            fvmop.template apply<stencil::StencilIloJlo>(lev);
            fvmop.template apply<stencil::StencilIloJhi>(lev);
        }
        if (!px && !pz) {
            fvmop.template apply<stencil::StencilIhiKlo>(lev);
            fvmop.template apply<stencil::StencilIhiKhi>(lev);
            fvmop.template apply<stencil::StencilIloKlo>(lev);
            fvmop.template apply<stencil::StencilIloKhi>(lev);
        }
        if (!py && !pz) {
            fvmop.template apply<stencil::StencilJhiKlo>(lev);
            fvmop.template apply<stencil::StencilJhiKhi>(lev);
            fvmop.template apply<stencil::StencilJloKlo>(lev);
            fvmop.template apply<stencil::StencilJloKhi>(lev);
        }

        // corners
        if (!px && !py && !pz) {
            fvmop.template apply<stencil::StencilIloJloKlo>(lev);
            fvmop.template apply<stencil::StencilIloJloKhi>(lev);
            fvmop.template apply<stencil::StencilIloJhiKlo>(lev);
            fvmop.template apply<stencil::StencilIloJhiKhi>(lev);
            fvmop.template apply<stencil::StencilIhiJloKlo>(lev);
            fvmop.template apply<stencil::StencilIhiJloKhi>(lev);
            fvmop.template apply<stencil::StencilIhiJhiKlo>(lev);
            fvmop.template apply<stencil::StencilIhiJhiKhi>(lev);
        }
    }
}

} // namespace kynema_sgf::fvm::impl

#endif /* FVM_UTILS_H */
