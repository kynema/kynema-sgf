#ifndef FVM_UTILS_H
#define FVM_UTILS_H

#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/fvm/stencils.H"

#include "AMReX_MFIter.H"

namespace amr_wind::fvm::impl {

/** Apply a finite volume operator for a given field
 */
template <typename FvmOp, typename FType>
AMREX_INLINE void apply(const FvmOp& fvmop, const FType& fld)
{
    namespace stencil = amr_wind::fvm::stencil;
    const int nlevels = fld.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& domain = fld.repo().mesh().Geom(lev).Domain();
        const auto& mfab = fld(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(mfab, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi) {

            fvmop.template apply<stencil::StencilInterior>(lev, mfi);

            // Check if the box touches any of the physical domain boundaries.
            // If not, short circuit the rest of the code and move on to the
            // next multifab.
            if (domain.strictly_contains(mfi.tilebox())) {
                continue;
            }

            // faces
            fvmop.template apply<stencil::StencilILO>(lev, mfi);
            fvmop.template apply<stencil::StencilJLO>(lev, mfi);
            fvmop.template apply<stencil::StencilKLO>(lev, mfi);
            fvmop.template apply<stencil::StencilIHI>(lev, mfi);
            fvmop.template apply<stencil::StencilJHI>(lev, mfi);
            fvmop.template apply<stencil::StencilKHI>(lev, mfi);

            // edges
            fvmop.template apply<stencil::StencilIhiJlo>(lev, mfi);
            fvmop.template apply<stencil::StencilIhiJhi>(lev, mfi);

            fvmop.template apply<stencil::StencilIhiKlo>(lev, mfi);
            fvmop.template apply<stencil::StencilIhiKhi>(lev, mfi);

            fvmop.template apply<stencil::StencilJhiKlo>(lev, mfi);
            fvmop.template apply<stencil::StencilJhiKhi>(lev, mfi);

            fvmop.template apply<stencil::StencilIloJlo>(lev, mfi);
            fvmop.template apply<stencil::StencilIloJhi>(lev, mfi);

            fvmop.template apply<stencil::StencilIloKlo>(lev, mfi);
            fvmop.template apply<stencil::StencilIloKhi>(lev, mfi);

            fvmop.template apply<stencil::StencilJloKlo>(lev, mfi);
            fvmop.template apply<stencil::StencilJloKhi>(lev, mfi);

            // corners
            fvmop.template apply<stencil::StencilIloJloKlo>(lev, mfi);
            fvmop.template apply<stencil::StencilIloJloKhi>(lev, mfi);
            fvmop.template apply<stencil::StencilIloJhiKlo>(lev, mfi);
            fvmop.template apply<stencil::StencilIloJhiKhi>(lev, mfi);
            fvmop.template apply<stencil::StencilIhiJloKlo>(lev, mfi);
            fvmop.template apply<stencil::StencilIhiJloKhi>(lev, mfi);
            fvmop.template apply<stencil::StencilIhiJhiKlo>(lev, mfi);
            fvmop.template apply<stencil::StencilIhiJhiKhi>(lev, mfi);
        }
    }
}

} // namespace amr_wind::fvm::impl

#endif /* FVM_UTILS_H */
