#ifndef HYDROSTATIC_OPS_H
#define HYDROSTATIC_OPS_H
#include <AMReX_MultiFabUtil.H>
#include "amr-wind/core/FieldRepo.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::hydrostatic {

static void define_rho0(
    amr_wind::Field& rho0,
    const amrex::Real rho1,
    const amrex::Real rho2,
    const amrex::Real wlev,
    const amrex::Vector<amrex::Geometry>& geom)
{
    for (int lev = 0; lev < rho0.repo().num_active_levels(); ++lev) {
        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();
        const auto& rho0_arrs = rho0(lev).arrays();
        amrex::ParallelFor(
            rho0(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real zbtm = problo[2] + (k * dx[2]);
                const amrex::Real vof = amrex::max<amrex::Real>(
                    amrex::min<amrex::Real>(1.0_rt, (wlev - zbtm) / dx[2]),
                    0.0_rt);
                rho0_arrs[nbx](i, j, k) =
                    (vof * rho1) + ((1.0_rt - vof) * rho2);
            });
    }
    amrex::Gpu::streamSynchronize();
}

static void define_p0(
    amr_wind::Field& p0,
    const amrex::Real rho1,
    const amrex::Real rho2,
    const amrex::Real wlev,
    const amrex::Real grav_z,
    const amrex::Vector<amrex::Geometry>& geom)
{
    for (int lev = 0; lev < p0.repo().num_active_levels(); ++lev) {
        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();
        const auto& probhi = geom[lev].ProbHiArray();
        auto p0_arrs = p0(lev).arrays();
        amrex::ParallelFor(
            p0(lev), p0.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real z = problo[2] + (k * dx[2]);
                // Integrated (top-down in z) phase heights to pressure node
                const amrex::Real ih_g = amrex::max<amrex::Real>(
                    0.0_rt, amrex::min(probhi[2] - wlev, probhi[2] - z));
                const amrex::Real ih_l = amrex::max<amrex::Real>(
                    0.0_rt,
                    amrex::min<amrex::Real>(wlev - z, wlev - problo[2]));
                // Integrated rho at pressure node
                const amrex::Real irho = (rho1 * ih_l) + (rho2 * ih_g);

                // Add term to reference pressure
                p0_arrs[nbx](i, j, k) = -irho * grav_z;
            });
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace amr_wind::hydrostatic

#endif
