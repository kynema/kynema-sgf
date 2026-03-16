#include "amr-wind/physics/multiphase/SloshingTank.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

SloshingTank::SloshingTank(CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
    , m_pressure(sim.repo().get_field("p"))
{
    amrex::ParmParse pp(identifier());
    pp.query("amplitude", m_amplitude);
    pp.query("peak_enhance", m_kappa);
    pp.query("water_level", m_waterlevel);

    pp.query("initialize_pressure", m_init_p);

    const auto& mphase = sim.physics_manager().get<MultiPhase>();
    m_rho1 = mphase.rho1();
    m_rho2 = mphase.rho2();
    m_gravity = mphase.gravity();
}

/** Initialize the velocity and levelset fields at the beginning of the
 *  simulation.
 */
void SloshingTank::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    velocity.setVal(0.0_rt);

    auto& levelset = m_levelset(level);
    auto& pressure = m_pressure(level);
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();
    const amrex::Real Amp = m_amplitude;
    const amrex::Real kappa = m_kappa;
    const amrex::Real water_level = m_waterlevel;
    const amrex::Real Lx = probhi[0] - problo[0];
    const amrex::Real Ly = probhi[1] - problo[1];
    const amrex::Real rho1 = m_rho1;
    const amrex::Real rho2 = m_rho2;
    const amrex::Real grav_z = m_gravity[2];

    const auto& phi_arrs = levelset.arrays();
    const auto& p = pressure.arrays();

    amrex::ParallelFor(
        levelset, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
            const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
            const amrex::Real z0 =
                water_level +
                (Amp * std::exp(
                           -kappa *
                           (std::pow(x - problo[0] - (0.5_rt * Lx), 2.0_rt) +
                            std::pow(y - problo[1] - (0.5_rt * Ly), 2.0_rt))));
            phi_arrs[nbx](i, j, k) = z0 - z;
        });

    if (m_init_p) {
        amrex::ParallelFor(
            pressure, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                // For pressure nodes, no offset
                const amrex::Real x = problo[0] + (i * dx[0]);
                const amrex::Real y = problo[1] + (j * dx[1]);
                const amrex::Real z = problo[2] + (k * dx[2]);
                const amrex::Real z0 =
                    water_level +
                    (Amp *
                     std::exp(
                         -kappa *
                         (std::pow(x - problo[0] - (0.5_rt * Lx), 2.0_rt) +
                          std::pow(y - problo[1] - (0.5_rt * Ly), 2.0_rt))));
                // Integrated (top-down in z) phase heights to pressure node
                amrex::Real ih_g = amrex::max<amrex::Real>(
                    0.0_rt,
                    amrex::min<amrex::Real>(probhi[2] - z0, probhi[2] - z));
                amrex::Real ih_l = amrex::max<amrex::Real>(
                    0.0_rt, amrex::min<amrex::Real>(z0 - z, z0 - problo[2]));
                // Integrated rho at pressure node
                const amrex::Real irho = (rho1 * ih_l) + (rho2 * ih_g);

                // Add term to reference pressure
                p[nbx](i, j, k) = -irho * grav_z;
            });
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace amr_wind
