#ifndef SPHERE_OPS_H
#define SPHERE_OPS_H

#include "amr-wind/immersed_boundary/bluff_body/Sphere.H"
#include "amr-wind/immersed_boundary/IBOps.H"
#include "amr-wind/immersed_boundary/IB.H"
#include "amr-wind/immersed_boundary/bluff_body/bluff_body_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::ib::ops {

template <>
struct ReadInputsOp<Sphere>
{
    void
    operator()(Sphere::DataType& data, const ::amr_wind::utils::MultiParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();

        bluff_body::read_inputs(wdata, info, pp);

        pp.get("center", wdata.center_loc);
        pp.get("radius", wdata.radius);

        amrex::Real search_radius = 3.0_rt * wdata.radius;
        // clang-format off
       const auto& origin=wdata.center_loc;
        info.bound_box = amrex::RealBox(
                origin[0] - search_radius,
                origin[1] - search_radius,
                origin[2] - search_radius,
                origin[0] + search_radius,
                origin[1] + search_radius,
                origin[2] + search_radius);
        // clang-format on
    }
};

template <>
struct InitDataOp<Sphere>
{
    void operator()(Sphere::DataType& data)
    {
        const auto& wdata = data.meta();

        auto& sim = data.sim();
        // cppcheck-suppress constVariableReference
        auto& mask_node = sim.repo().get_int_field("mask_node");
        // cppcheck-suppress constVariableReference
        auto& levelset = sim.repo().get_field("ib_levelset");

        auto nlevels = sim.repo().num_active_levels();
        auto geom = sim.mesh().Geom();

        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& problo = geom[lev].ProbLoArray();
            const auto& dx = geom[lev].CellSizeArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for (amrex::MFIter mfi(levelset(lev)); mfi.isValid(); ++mfi) {
                const auto& bx = mfi.growntilebox();
                const auto& phi = levelset(lev).array(mfi);

                const amrex::Real x0 = wdata.center_loc[0];
                const amrex::Real y0 = wdata.center_loc[1];
                const amrex::Real z0 = wdata.center_loc[2];
                const amrex::Real R = wdata.radius;
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        const amrex::Real x =
                            problo[0] + ((i + 0.5_rt) * dx[0]);
                        const amrex::Real y =
                            problo[1] + ((j + 0.5_rt) * dx[1]);
                        const amrex::Real z =
                            problo[2] + ((k + 0.5_rt) * dx[2]);
                        const amrex::Real phi_glob = phi(i, j, k);
                        const amrex::Real r = std::sqrt(
                            ((x - x0) * (x - x0)) + ((y - y0) * (y - y0)) +
                            ((z - z0) * (z - z0)));

                        const amrex::Real phi_loc = r - R;
                        phi(i, j, k) =
                            amrex::min<amrex::Real>(phi_loc, phi_glob);
                    });

                const auto& nbx = mfi.nodaltilebox();
                const auto& epsilon_node = mask_node(lev).array(mfi);
                amrex::ParallelFor(
                    nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        const amrex::Real x = problo[0] + (i * dx[0]);
                        const amrex::Real y = problo[1] + (j * dx[1]);
                        const amrex::Real z = problo[2] + (k * dx[2]);

                        const amrex::Real r = std::sqrt(
                            ((x - x0) * (x - x0)) + ((y - y0) * (y - y0)) +
                            ((z - z0) * (z - z0)));
                        if (r <= R) {
                            epsilon_node(i, j, k) = 0;
                        }
                    });
            }
        }
    }
};

} // namespace amr_wind::ib::ops

#endif /* CYLINDER_OPS_H */
