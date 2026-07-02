#include "src/physics/multiphase/FreeSurfaceDamping.H"
#include "src/CFDSim.H"
#include "src/core/FieldRepo.H"
#include "src/core/MultiParser.H"
#include "src/physics/multiphase/MultiPhase.H"
#include "src/utilities/IOManager.H"
#include "src/ocean_waves/utils/wave_utils_K.H"

#include <algorithm>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf {

FreeSurfaceDamping::FreeSurfaceDamping(CFDSim& sim)
    : m_repo(sim.repo()), m_time(sim.time())
{
    amrex::ParmParse pp(identifier());

    pp.query("vertical_velocity_damping", m_vertical_velocity_damping);
    pp.query("volume_fraction_damping", m_volume_fraction_damping);

    pp.query("global_damping", m_global_damping);

    pp.query("time_scale_fraction", m_time_scale);
    if (m_time_scale < 0.0_rt || m_time_scale > 1.0_rt) {
        amrex::Abort(
            "FreeSurfaceDamping.time_scale_fraction must be in [0, 1].");
    }

    if (!m_global_damping) {
        pp.get("length_xlo", m_length_xlo);
        pp.get("length_xhi", m_length_xhi);
        pp.get("length_ylo", m_length_ylo);
        pp.get("length_yhi", m_length_yhi);
    }

    const auto& mphase = sim.physics_manager().get<MultiPhase>();
    m_rho1 = mphase.rho1();
    m_rho2 = mphase.rho2();
}

FreeSurfaceDamping::~FreeSurfaceDamping() = default;

void FreeSurfaceDamping::post_advance_work()
{
    BL_PROFILE("kynema-sgf::FreeSurfaceDamping::post_advance_work");

    const int nlevels = m_repo.num_active_levels();
    auto& vof = m_repo.get_field("vof");
    auto& velocity = m_repo.get_field("velocity");
    auto& density = m_repo.get_field("density");

    constexpr amrex::Real vof_tiny = 1.0e-12_rt;
    const auto vertical_velocity_damping = m_vertical_velocity_damping;
    const auto volume_fraction_damping = m_volume_fraction_damping;
    const auto global_damping = m_global_damping;
    const auto l_xlo = m_length_xlo;
    const auto l_xhi = m_length_xhi;
    const auto l_ylo = m_length_ylo;
    const auto l_yhi = m_length_yhi;
    const auto time_scale = m_time_scale;
    const auto rho1 = m_rho1;
    const auto rho2 = m_rho2;

    // Get time
    const auto& time = m_time.new_time();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = m_repo.mesh().Geom(lev).CellSizeArray();
        const auto& problo = m_repo.mesh().Geom(lev).ProbLoArray();
        const auto& probhi = m_repo.mesh().Geom(lev).ProbHiArray();
        auto vel_arrs = velocity(lev).arrays();
        auto rho_arrs = density(lev).arrays();
        auto volfrac_arrs = vof(lev).arrays();

        amrex::ParallelFor(
            velocity(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);

                auto vel = vel_arrs[nbx];
                auto rho = rho_arrs[nbx];
                auto volfrac = volfrac_arrs[nbx];

                // Initialize Gamma as 0, turning damping it on globally
                Gamma = 0.0_rt;

                if (!global_damping) {
                    // If not global, get gamma for each possible direction
                    const amrex::Real Gamma_xlo =
                        l_xlo > constants::EPS
                            ? ocean_waves::utils::gamma_generate(
                                  x - problo[0], l_xlo)
                            : 1.0_rt;
                    const amrex::Real Gamma_xhi =
                        l_xhi > constants::EPS
                            ? ocean_waves::utils::gamma_absorb(
                                  x - (probhi[0] - l_xhi), l_xhi, 1.0_rt)
                            : 1.0_rt;
                    const amrex::Real Gamma_ylo =
                        l_ylo > constants::EPS
                            ? ocean_waves::utils::gamma_generate(
                                  y - problo[1], l_ylo)
                            : 1.0_rt;
                    const amrex::Real Gamma_yhi =
                        l_yhi > constants::EPS
                            ? ocean_waves::utils::gamma_absorb(
                                  y - (probhi[1] - l_yhi), l_yhi, 1.0_rt)
                            : 1.0_rt;

                    Gamma = amrex::min<amrex::Real>(
                        amrex::min<amrex::Real>(Gamma_xhi, Gamma_xlo),
                        amrex::min<amrex::Real>(Gamma_yhi, Gamma_ylo));
                }

                // Skip if Gamma is close enough to 1
                bool outside_zones = Gamma + constants::EPS >= 1.0_rt;

                // Determine if cell is multiphase
                bool is_multiphase =
                    (volfrac(i, j, k) > vof_tiny &&
                     volfrac(i, j, k) < 1.0_rt - vof_tiny);

                if (!outside_zones && is_multiphase) {

                    // Vertical velocity damping
                    // - Damps velocity toward zero
                    if (vertical_velocity_damping) {
                        const amrex::Real dvel =
                            ocean_waves::utils::combine_linear(
                                Gamma, 0.0_rt, vel(i, j, k, 2)) -
                            vel(i, j, k, 2);
                        vel(i, j, k, 2) += dvel * time_scale;
                    }

                    // Volume fraction damping
                    // - Smooths volume fraction based on horizontal average,
                    // avoiding an assumption of target water level
                    if (volume_fraction_damping) {
                        const amrex::Real smooth_vof =
                            (volfrac(i, j, k) + volfrac(i - 1, j, k) +
                             volfrac(i + 1, j, k) + volfrac(i, j - 1, k) +
                             volfrac(i, j + 1, k) + volfrac(i - 1, j - 1, k) +
                             volfrac(i + 1, j + 1, k) +
                             volfrac(i - 1, j + 1, k) +
                             volfrac(i + 1, j - 1, k)) /
                            9.0_rt;
                        amrex::Real dvof =
                            ocean_waves::utils::combine_linear(
                                Gamma, smooth_vof, volfrac(i, j, k)) -
                            volfrac(i, j, k);
                        volfrac(i, j, k) += dvof * time_scale;

                        if (amrex::Math::abs(dvof) > constants::EPS) {
                            // Update density based on change in volume fraction
                            rho(i, j, k) = (rho1 * volfrac(i, j, k)) +
                                           (rho2 * (1.0_rt - volfrac(i, j, k)));
                        }
                    }
                }
            });
    }
    amrex::Gpu::streamSynchronize();

    if (vertical_velocity_damping) {
        velocity.fillpatch(time);
    }
    if (volume_fraction_damping) {
        vof.fillpatch(time);
        density.fillpatch(time);
    }
}

} // namespace kynema_sgf
