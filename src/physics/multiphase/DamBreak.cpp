#include <numbers>
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/physics/multiphase/DamBreak.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

DamBreak::DamBreak(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
    , m_density(sim.repo().get_field("density"))
{
    amrex::ParmParse pp(identifier());
    pp.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp.query("width", m_width);
    pp.query("height", m_height);
}

/** Initialize the velocity and levelset fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::DamBreakFieldInit
 */
void DamBreak::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    velocity.setVal(0.0_rt, 0, AMREX_SPACEDIM);

    auto& levelset = m_levelset(level);
    auto& density = m_density(level);
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const amrex::Real xc = m_loc[0];
    const amrex::Real zc = m_loc[2];
    const amrex::Real width = m_width;
    const amrex::Real height = m_height;

    const auto& mphase = m_sim.physics_manager().get<MultiPhase>();
    const amrex::Real rho1 = mphase.rho1();
    const amrex::Real rho2 = mphase.rho2();

    const auto& phi_arrs = levelset.arrays();
    const auto& rho_arrs = density.arrays();
    const amrex::Real eps = std::cbrt(2.0_rt * dx[0] * dx[1] * dx[2]);

    amrex::ParallelFor(
        levelset, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

            if (x - xc < width && z - zc < height) {
                phi_arrs[nbx](i, j, k) = amrex::min(width - x, height - z);
            } else if (x - xc < width && z - zc > height) {
                phi_arrs[nbx](i, j, k) = height - (z - zc);
            } else if (x - xc > width && z - zc < height) {
                phi_arrs[nbx](i, j, k) = width - (x - xc);
            } else {
                phi_arrs[nbx](i, j, k) =
                    amrex::min(width - (x - xc), height - (z - zc));
            }
            amrex::Real smooth_heaviside;
            if (phi_arrs[nbx](i, j, k) > eps) {
                smooth_heaviside = 1.0_rt;
            } else if (phi_arrs[nbx](i, j, k) < -eps) {
                smooth_heaviside = 0.0_rt;
            } else {
                smooth_heaviside =
                    0.5_rt * (1.0_rt + phi_arrs[nbx](i, j, k) / eps +
                              1.0_rt / std::numbers::pi_v<amrex::Real> *
                                  std::sin(
                                      phi_arrs[nbx](i, j, k) *
                                      std::numbers::pi_v<amrex::Real> / eps));
            }
            rho_arrs[nbx](i, j, k) = (rho1 * smooth_heaviside) +
                                     (rho2 * (1.0_rt - smooth_heaviside));
        });
    amrex::Gpu::streamSynchronize();
}

} // namespace amr_wind
