#include <numbers>
#include "amr-wind/immersed_boundary/bluff_body/bluff_body_ops.H"
#include "amr-wind/core/MultiParser.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"
#include "amr-wind/utilities/io_utils.H"

#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"

// Used for mms
#include "amr-wind/physics/ConvectingTaylorVortex.H"

#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::ib::bluff_body {

void read_inputs(
    BluffBodyBaseData& wdata,
    IBInfo& /*unused*/,
    const ::amr_wind::utils::MultiParser& pp)
{
    pp.query("has_wall_model", wdata.has_wall_model);
    pp.query("is_moving", wdata.is_moving);
    pp.query("is_mms", wdata.is_mms);
    pp.queryarr("vel_bc", wdata.vel_bc);
}

void init_data_structures(BluffBodyBaseData& /*unused*/) {}

void apply_mms_vel(CFDSim& sim)
{
    const int nlevels = sim.repo().num_active_levels();

    const auto& levelset = sim.repo().get_field("ib_levelset");
    auto& velocity = sim.repo().get_field("velocity");
    auto& m_conv_taylor_green =
        sim.physics_manager().get<ctv::ConvectingTaylorVortex>();

    const amrex::Real u0 = m_conv_taylor_green.get_u0();
    const amrex::Real v0 = m_conv_taylor_green.get_v0();
    const amrex::Real omega = m_conv_taylor_green.get_omega();

    amrex::Real t = sim.time().new_time();
    auto& geom = sim.mesh().Geom();

    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();

        const auto& phi_arrs = levelset(lev).const_arrays();
        const auto& varrs = velocity(lev).arrays();
        amrex::ParallelFor(
            levelset(lev), levelset.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);

                if (phi_arrs[nbx](i, j, k) <= 0) {
                    varrs[nbx](i, j, k, 0) =
                        u0 -
                        (std::cos(
                             std::numbers::pi_v<amrex::Real> * (x - u0 * t)) *
                         std::sin(
                             std::numbers::pi_v<amrex::Real> * (y - v0 * t)) *
                         std::exp(-2.0_rt * omega * t));
                    varrs[nbx](i, j, k, 1) =
                        v0 +
                        (std::sin(
                             std::numbers::pi_v<amrex::Real> * (x - u0 * t)) *
                         std::cos(
                             std::numbers::pi_v<amrex::Real> * (y - v0 * t)) *
                         std::exp(-2.0_rt * omega * t));
                    varrs[nbx](i, j, k, 2) = 0.0_rt;
                }
            });
    }
    amrex::Gpu::streamSynchronize();
}

void apply_dirichlet_vel(CFDSim& sim, const amrex::Vector<amrex::Real>& vel_bc)
{
    const int nlevels = sim.repo().num_active_levels();
    auto& geom = sim.mesh().Geom();
    auto& velocity = sim.repo().get_field("velocity");
    auto& levelset = sim.repo().get_field("ib_levelset");
    levelset.fillpatch(sim.time().current_time());
    auto& normal = sim.repo().get_field("ib_normal");
    fvm::gradient(normal, levelset);
    field_ops::normalize(normal);
    normal.fillpatch(sim.time().current_time());

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = geom[lev].CellSizeArray();
        // const auto& problo = geom[lev].ProbLoArray();
        // Defining the "ghost-cell" band distance
        amrex::Real phi_b = std::cbrt(dx[0] * dx[1] * dx[2]);

        const auto& varrs = velocity(lev).arrays();
        const auto& phi_arrs = levelset(lev).const_arrays();
        const auto& norm_arrs = normal(lev).arrays();

        const amrex::Real velx = vel_bc[0];
        const amrex::Real vely = vel_bc[1];
        const amrex::Real velz = vel_bc[2];

        amrex::ParallelFor(
            levelset(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                // Pure solid-body points
                if (phi_arrs[nbx](i, j, k) < -phi_b) {
                    varrs[nbx](i, j, k, 0) = velx;
                    varrs[nbx](i, j, k, 1) = vely;
                    varrs[nbx](i, j, k, 2) = velz;
                    norm_arrs[nbx](i, j, k, 0) = 0.;
                    norm_arrs[nbx](i, j, k, 1) = 0.;
                    norm_arrs[nbx](i, j, k, 2) = 0.;
                    // This determines the ghost-cells
                } else if (
                    phi_arrs[nbx](i, j, k) < 0 &&
                    phi_arrs[nbx](i, j, k) >= -phi_b) {
                    // For this particular ghost-cell find the
                    // body-intercept (BI) point and image-point (IP)
                    // First define the ghost cell point
                    // amrex::Real x_GC = problo[0] + (i + 0.5_rt) * dx[0];
                    // amrex::Real y_GC = problo[1] + (j + 0.5_rt) * dx[1];
                    // amrex::Real z_GC = problo[2] + (k + 0.5_rt) * dx[2];
                    // Find the "image-points"
                    varrs[nbx](i, j, k, 0) = velx;
                    varrs[nbx](i, j, k, 1) = vely;
                    varrs[nbx](i, j, k, 2) = velz;
                } else {
                    norm_arrs[nbx](i, j, k, 0) = 0.;
                    norm_arrs[nbx](i, j, k, 1) = 0.;
                    norm_arrs[nbx](i, j, k, 2) = 0.;
                }
            });
    }
    amrex::Gpu::streamSynchronize();
}

void prepare_netcdf_file(
    const std::string& ncfile,
    const BluffBodyBaseData& meta,
    const IBInfo& info)
{
    amrex::ignore_unused(ncfile, meta, info);
}

void write_netcdf(
    const std::string& ncfile,
    const BluffBodyBaseData& meta,
    const IBInfo& info,
    const amrex::Real time)
{
    amrex::ignore_unused(ncfile, meta, info, time);
}

} // namespace amr_wind::ib::bluff_body
