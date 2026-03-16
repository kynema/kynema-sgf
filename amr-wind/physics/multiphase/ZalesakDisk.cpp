#include <numbers>
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/physics/multiphase/ZalesakDisk.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/core/field_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

ZalesakDisk::ZalesakDisk(CFDSim& sim)
    : m_sim(sim)
    , m_velocity(sim.repo().get_field("velocity"))
    , m_levelset(sim.repo().get_field("levelset"))
    , m_density(sim.repo().get_field("density"))
{
    amrex::ParmParse pp(identifier());
    pp.queryarr("location", m_loc, 0, AMREX_SPACEDIM);
    pp.query("radius", m_radius);
    pp.query("period", m_TT);
    amrex::ParmParse pinc("incflo");
    pinc.add("prescribe_velocity", true);
}

/** Initialize the velocity and levelset fields at the beginning of the
 *  simulation.
 *
 *  \sa amr_wind::ZalesakDiskFieldInit
 */
void ZalesakDisk::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& levelset = m_levelset(level);
    auto& density = m_density(level);
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();

    const auto& mphase = m_sim.physics_manager().get<MultiPhase>();
    const amrex::Real rho1 = mphase.rho1();
    const amrex::Real rho2 = mphase.rho2();

    auto& u_mac = m_sim.repo().get_field("u_mac")(level);
    auto& v_mac = m_sim.repo().get_field("v_mac")(level);
    auto& w_mac = m_sim.repo().get_field("w_mac")(level);

    const amrex::Real xc = m_loc[0];
    const amrex::Real yc = m_loc[1];
    const amrex::Real zc = m_loc[2];
    const amrex::Real radius = m_radius;
    const amrex::Real TT = m_TT;
    const amrex::Real hwidth = m_halfwidth;
    const amrex::Real depth = m_depth;
    const auto& uf_arrs = u_mac.arrays();
    const auto& vf_arrs = v_mac.arrays();
    const auto& wf_arrs = w_mac.arrays();
    const auto& vel_arrs = velocity.arrays();
    const auto& phi_arrs = levelset.arrays();
    const auto& rho_arrs = density.arrays();

    amrex::ParallelFor(
        levelset, amrex::IntVect(1),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
            const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

            uf_arrs[nbx](i, j, k) =
                2.0_rt * std::numbers::pi_v<amrex::Real> / TT * (0.5_rt - y);
            vf_arrs[nbx](i, j, k) =
                2.0_rt * std::numbers::pi_v<amrex::Real> / TT * (x - 0.5_rt);
            wf_arrs[nbx](i, j, k) = 0.0_rt;

            vel_arrs[nbx](i, j, k, 0) =
                2.0_rt * std::numbers::pi_v<amrex::Real> / TT * (0.5_rt - y);
            vel_arrs[nbx](i, j, k, 1) =
                2.0_rt * std::numbers::pi_v<amrex::Real> / TT * (x - 0.5_rt);
            vel_arrs[nbx](i, j, k, 2) = 0.0_rt;

            // First define the sphere
            const amrex::Real r = std::sqrt(
                ((x - xc) * (x - xc)) + ((y - yc) * (y - yc)) +
                ((z - zc) * (z - zc)));
            phi_arrs[nbx](i, j, k) = radius - r;

            // Then the slot
            // Signed distances in lateral (x, y) directions
            const amrex::Real sd_xr = -hwidth + (x - xc);
            const amrex::Real sd_xl = -hwidth - (x - xc);
            const amrex::Real sd_x = amrex::max<amrex::Real>(sd_xr, sd_xl);

            const amrex::Real sd_y = radius - depth - (y - yc);
            const amrex::Real min_signed_dist =
                amrex::max<amrex::Real>(sd_x, sd_y);

            // Additional distance if past sphere (distance to corners)
            const amrex::Real reduced_radius =
                std::sqrt((radius * radius) - (hwidth * hwidth));
            const amrex::Real r_2D =
                std::sqrt(std::pow(y - yc, 2.0_rt) + std::pow(z - zc, 2.0_rt));
            const amrex::Real sd_r = -std::sqrt(
                std::pow(r_2D - reduced_radius, 2.0_rt) +
                std::pow(sd_x, 2.0_rt));

            const bool in_slot_x_ymin =
                y - yc > radius - depth && std::abs(x - xc) < hwidth;
            const bool in_slot_r = r_2D < reduced_radius;

            if (in_slot_x_ymin) {
                // Prescribe slot distances directly (overwrite sphere)
                if (in_slot_r) {
                    phi_arrs[nbx](i, j, k) = min_signed_dist;
                } else {
                    phi_arrs[nbx](i, j, k) = sd_r;
                }
            } else {
                // Select the minimum of the two
                phi_arrs[nbx](i, j, k) =
                    amrex::min(phi_arrs[nbx](i, j, k), min_signed_dist);
            }

            amrex::Real smooth_heaviside;
            const amrex::Real eps = std::cbrt(2.0_rt * dx[0] * dx[1] * dx[2]);
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

    m_levelset.fillpatch(m_sim.time().current_time());
    m_velocity.fillpatch(m_sim.time().current_time());
    m_density.fillpatch(m_sim.time().current_time());
}

void ZalesakDisk::pre_advance_work()
{

    const int nlevels = m_sim.repo().num_active_levels();
    const auto& geom = m_sim.mesh().Geom();

    // Overriding the velocity field
    for (int lev = 0; lev < nlevels; ++lev) {
        auto& u_mac = m_sim.repo().get_field("u_mac")(lev);
        auto& v_mac = m_sim.repo().get_field("v_mac")(lev);
        auto& w_mac = m_sim.repo().get_field("w_mac")(lev);
        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();
        const amrex::Real TT = m_TT;
        const auto& uf_arrs = u_mac.arrays();
        const auto& vf_arrs = v_mac.arrays();
        const auto& wf_arrs = w_mac.arrays();
        amrex::ParallelFor(
            m_velocity(lev), amrex::IntVect(1),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);

                uf_arrs[nbx](i, j, k) = 2.0_rt *
                                        std::numbers::pi_v<amrex::Real> / TT *
                                        (0.5_rt - y);
                vf_arrs[nbx](i, j, k) = 2.0_rt *
                                        std::numbers::pi_v<amrex::Real> / TT *
                                        (x - 0.5_rt);
                wf_arrs[nbx](i, j, k) = 0.0_rt;
            });
        amrex::Gpu::streamSynchronize();
        u_mac.FillBoundary(geom[lev].periodicity());
        v_mac.FillBoundary(geom[lev].periodicity());
        w_mac.FillBoundary(geom[lev].periodicity());
    }
}

void ZalesakDisk::post_advance_work()
{

    const int nlevels = m_sim.repo().num_active_levels();
    const auto& geom = m_sim.mesh().Geom();

    // Overriding the velocity field
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();
        const amrex::Real TT = m_TT;
        const auto& vel_arrs = m_velocity(lev).arrays();
        amrex::ParallelFor(
            m_velocity(lev),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);

                vel_arrs[nbx](i, j, k, 0) = 2.0_rt *
                                            std::numbers::pi_v<amrex::Real> /
                                            TT * (0.5_rt - y);
                vel_arrs[nbx](i, j, k, 1) = 2.0_rt *
                                            std::numbers::pi_v<amrex::Real> /
                                            TT * (x - 0.5_rt);
                vel_arrs[nbx](i, j, k, 2) = 0.0_rt;
            });
    }
    amrex::Gpu::streamSynchronize();

    m_velocity.fillpatch(m_sim.time().current_time());
}

} // namespace amr_wind
