#include <numbers>
#include "amr-wind/physics/multiphase/RainDrop.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

RainDrop::RainDrop(CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
    , m_sim(sim)
{
    amrex::ParmParse pp(identifier());
    pp.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp.query("radius", m_radius);
    pp.queryarr("droplet_velocity", m_vel, 0, AMREX_SPACEDIM);
}

/** Initialize the velocity and levelset fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::RainDropFieldInit
 */
void RainDrop::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& levelset = m_levelset(level);
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();

    const amrex::Real xc = m_loc[0];
    const amrex::Real yc = m_loc[1];
    const amrex::Real zc = m_loc[2];
    const amrex::Real radius = m_radius;
    const auto& phi_arrs = levelset.arrays();
    amrex::ParallelFor(
        velocity, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
            const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

            // Calculate level-set
            phi_arrs[nbx](i, j, k) =
                radius - std::sqrt(
                             ((x - xc) * (x - xc)) + ((y - yc) * (y - yc)) +
                             ((z - zc) * (z - zc)));
        });
    amrex::Gpu::streamSynchronize();

    auto& mphase = m_sim.physics_manager().get<MultiPhase>();
    // Get VOF array if it is valid
    const bool yes_vof = m_sim.repo().field_exists("vof");
    amrex::MultiFab* work_ptr;
    if (yes_vof) {
        auto& volfrac = m_sim.repo().get_field("vof")(level);
        mphase.levelset2vof();
        work_ptr = &volfrac;
    } else {
        work_ptr = &levelset;
    }

    // Use density to calculate mixed velocity in multiphase cells
    const amrex::Real rhol = mphase.rho1();
    const amrex::Real rhog = mphase.rho2();
    const amrex::Real eps = std::cbrt(2.0_rt * dx[0] * dx[1] * dx[2]);
    amrex::Gpu::DeviceVector<amrex::Real> dvel(AMREX_SPACEDIM);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, m_vel.begin(), m_vel.end(), dvel.begin());
    const auto* vptr = dvel.data();
    const auto& vel_arrs = velocity.arrays();
    const auto& work_arrs = work_ptr->arrays();
    amrex::ParallelFor(
        velocity, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            // Calculate mass-weighted velocity (gas vel is 0)
            amrex::Real vof;
            if (yes_vof) {
                // Work array here is vof
                vof = work_arrs[nbx](i, j, k);
            } else {
                // Work array here is phi
                if (work_arrs[nbx](i, j, k) > eps) {
                    vof = 1.0_rt;
                } else if (work_arrs[nbx](i, j, k) < -eps) {
                    vof = 0.0_rt;
                } else {
                    vof = 0.5_rt *
                          (1.0_rt + work_arrs[nbx](i, j, k) / eps +
                           1.0_rt / std::numbers::pi_v<amrex::Real> *
                               std::sin(
                                   work_arrs[nbx](i, j, k) *
                                   std::numbers::pi_v<amrex::Real> / eps));
                }
            }
            const amrex::Real dens = ((1.0_rt - vof) * rhog) + (vof * rhol);
            vel_arrs[nbx](i, j, k, 0) = vof * rhol * vptr[0] / dens;
            vel_arrs[nbx](i, j, k, 1) = vof * rhol * vptr[1] / dens;
            vel_arrs[nbx](i, j, k, 2) = vof * rhol * vptr[2] / dens;
        });

    amrex::Gpu::streamSynchronize();
}

} // namespace amr_wind
