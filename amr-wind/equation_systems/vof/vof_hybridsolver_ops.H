#ifndef VOF_HYBRIDSOLVER_OPS_H
#define VOF_HYBRIDSOLVER_OPS_H

#include <cmath>
#include "amr-wind/core/FieldRepo.H"
#include "AMReX_FArrayBox.H"

namespace amr_wind::multiphase {

static void replace_masked_vof(
    const int nlevels,
    amr_wind::IntField& f_iblank,
    amr_wind::Field& f_vof,
    amr_wind::Field& f_vof_new)
{
    // Sharpen data from nalu-wind (in iblank regions)
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& native_flag_arrs = f_iblank(lev).const_arrays();
        const auto& volfrac_arrs = f_vof(lev).arrays();
        const auto& vfmasked_arrs = f_vof_new(lev).const_arrays();
        amrex::ParallelFor(
            f_iblank(lev), f_iblank.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                // In iblanked regions, sharpen VOF and limit it
                volfrac_arrs[nbx](i, j, k) =
                    (native_flag_arrs[nbx](i, j, k) > 0)
                        ? volfrac_arrs[nbx](i, j, k)
                        : vfmasked_arrs[nbx](i, j, k);
            });
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace amr_wind::multiphase

#endif // VOF_HYBRIDSOLVER_OPS.H
