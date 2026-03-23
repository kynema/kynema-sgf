#include <numbers>
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/physics/multiphase/BreakingWaves.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

BreakingWaves::BreakingWaves(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
    , m_density(sim.repo().get_field("density"))
{

    {
        std::string turbulence_model;
        amrex::ParmParse pp("turbulence");
        pp.query("model", turbulence_model);
        if (turbulence_model == "Laminar") {
            m_laminar = true;
        }
    }

    {
        amrex::ParmParse pp("BreakingWaves");
        pp.query("amplitude", m_amplitude);
        pp.query("wavelength", m_wavelength);
        pp.query("water_surface_level", m_waterlevel);

        if (m_laminar) {
            pp.query("air_vel_shear_rate", m_air_vel_shear_rate);
        } else {
            pp.query("air_vel_shear_rate", m_air_vel_shear_rate);
            pp.query("tke0", m_tke0);
            pp.query("sdr0", m_sdr0);
        }
    }
}

BreakingWaves::~BreakingWaves() = default;

/** Initialize the velocity and vof fields at the beginning of the
 *  simulation.
 */
void BreakingWaves::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& levelset = m_levelset(level);
    auto& density = m_density(level);

    const auto& mphase = m_sim.physics_manager().get<MultiPhase>();
    const amrex::Real rho1 = mphase.rho1();
    const amrex::Real rho2 = mphase.rho2();
    velocity.setVal(0.0_rt, 0, AMREX_SPACEDIM);

    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const amrex::Real alpha = m_amplitude;
    const amrex::Real lambda = m_wavelength;
    const amrex::Real water_level = m_waterlevel;
    const amrex::Real vel_shear_rate = m_air_vel_shear_rate;

    if (!m_laminar) {
        auto& tke = m_sim.repo().get_field("tke")(level);
        auto& sdr = m_sim.repo().get_field("sdr")(level);
        tke.setVal(m_tke0);
        sdr.setVal(m_sdr0);
    }

    const auto& vel_arrs = velocity.arrays();
    const auto& phi_arrs = levelset.arrays();
    const auto& rho_arrs = density.arrays();

    amrex::ParallelFor(
        levelset, m_levelset.num_grow(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
            const amrex::Real kappa =
                2.0_rt * std::numbers::pi_v<amrex::Real> / lambda;
            const amrex::Real epsilon = alpha * kappa;
            // Compute free surface amplitude
            const amrex::Real eta =
                water_level +
                (alpha * ((1.0_rt - 1.0_rt / 16.0_rt * epsilon * epsilon) *
                              std::cos(kappa * x) +
                          0.5_rt * epsilon * std::cos(2.0_rt * kappa * x) +
                          3.0_rt / 8.0_rt * epsilon * epsilon *
                              std::cos(3.0_rt * kappa * x)));
            phi_arrs[nbx](i, j, k) = eta - z;
            // compute velocities
            const amrex::Real g = 9.81_rt;
            const amrex::Real Omega =
                std::sqrt(g * kappa * (1.0_rt + epsilon * epsilon));
            if (z < eta) {
                vel_arrs[nbx](i, j, k, 0) =
                    Omega * alpha * std::exp(kappa * z) * std::cos(kappa * x);
                vel_arrs[nbx](i, j, k, 2) =
                    Omega * alpha * std::exp(kappa * z) * std::sin(kappa * x);
            } else {
                vel_arrs[nbx](i, j, k, 0) = vel_shear_rate * (z - eta);
            }
            // compute density
            amrex::Real smooth_heaviside;
            amrex::Real eps = std::cbrt(2.0_rt * dx[0] * dx[1] * dx[2]);
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
