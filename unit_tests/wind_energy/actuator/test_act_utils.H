#ifndef TEST_ACT_UTILS_H
#define TEST_ACT_UTILS_H

#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

// Utility function to populate the velocity field used for tests
AMREX_INLINE void init_field(amr_wind::Field& fld)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();
    const int ncomp = fld.num_comp();

    amrex::Real offsetx = 0.0_rt;
    amrex::Real offsety = 0.0_rt;
    amrex::Real offsetz = 0.0_rt;
    if (fld.field_location() == amr_wind::FieldLoc::CELL) {
        offsetx = 0.5_rt;
        offsety = 0.5_rt;
        offsetz = 0.5_rt;
    }
    if (fld.field_location() == amr_wind::FieldLoc::XFACE) {
        offsety = 0.5_rt;
        offsetz = 0.5_rt;
    }
    if (fld.field_location() == amr_wind::FieldLoc::YFACE) {
        offsetx = 0.5_rt;
        offsetz = 0.5_rt;
    }
    if (fld.field_location() == amr_wind::FieldLoc::ZFACE) {
        offsetx = 0.5_rt;
        offsety = 0.5_rt;
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();

        const auto& farrs = fld(lev).arrays();
        amrex::ParallelFor(
            fld(lev), fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real x = problo[0] + ((i + offsetx) * dx[0]);
                const amrex::Real y = problo[1] + ((j + offsety) * dx[1]);
                const amrex::Real z = problo[2] + ((k + offsetz) * dx[2]);

                for (int d = 0; d < ncomp; d++) {
                    farrs[nbx](i, j, k, d) = x + y + z;
                }
            });
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace amr_wind_tests

#endif /* TEST_ACT_UTILS_H */
