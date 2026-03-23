#include <numbers>
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/physics/multiphase/VortexPatch.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

VortexPatch::VortexPatch(CFDSim& sim)
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
 *  \sa amr_wind::VortexPatchFieldInit
 */
void VortexPatch::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& levelset = m_levelset(level);
    auto& density = m_density(level);
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();

    const amrex::Real xc = m_loc[0];
    const amrex::Real yc = m_loc[1];
    const amrex::Real zc = m_loc[2];
    const amrex::Real radius = m_radius;

    const auto& mphase = m_sim.physics_manager().get<MultiPhase>();
    const amrex::Real rho1 = mphase.rho1();
    const amrex::Real rho2 = mphase.rho2();

    auto& u_mac = m_sim.repo().get_field("u_mac")(level);
    auto& v_mac = m_sim.repo().get_field("v_mac")(level);
    auto& w_mac = m_sim.repo().get_field("w_mac")(level);

    const auto& uf_arrs = u_mac.arrays();
    const auto& vf_arrs = v_mac.arrays();
    const auto& wf_arrs = w_mac.arrays();
    const auto& vel_arrs = velocity.arrays();
    const auto& phi_arrs = levelset.arrays();
    const auto& rho_arrs = density.arrays();
    const amrex::Real eps = std::cbrt(2.0_rt * dx[0] * dx[1] * dx[2]);

    amrex::ParallelFor(
        velocity, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
            const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
            const amrex::Real xf = problo[0] + (i * dx[0]);
            const amrex::Real yf = problo[1] + (j * dx[1]);
            const amrex::Real zf = problo[2] + (k * dx[2]);
            uf_arrs[nbx](i, j, k) =
                2.0_rt * std::sin(std::numbers::pi_v<amrex::Real> * xf) *
                std::sin(std::numbers::pi_v<amrex::Real> * xf) *
                std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * y) *
                std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * z);
            vf_arrs[nbx](i, j, k) =
                -std::sin(std::numbers::pi_v<amrex::Real> * yf) *
                std::sin(std::numbers::pi_v<amrex::Real> * yf) *
                std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * x) *
                std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * z);
            wf_arrs[nbx](i, j, k) =
                -std::sin(std::numbers::pi_v<amrex::Real> * zf) *
                std::sin(std::numbers::pi_v<amrex::Real> * zf) *
                std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * x) *
                std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * y);

            vel_arrs[nbx](i, j, k, 0) =
                2.0_rt * std::sin(std::numbers::pi_v<amrex::Real> * x) *
                std::sin(std::numbers::pi_v<amrex::Real> * x) *
                std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * y) *
                std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * z);
            vel_arrs[nbx](i, j, k, 1) =
                -std::sin(std::numbers::pi_v<amrex::Real> * y) *
                std::sin(std::numbers::pi_v<amrex::Real> * y) *
                std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * x) *
                std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * z);
            vel_arrs[nbx](i, j, k, 2) =
                -std::sin(std::numbers::pi_v<amrex::Real> * z) *
                std::sin(std::numbers::pi_v<amrex::Real> * z) *
                std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * x) *
                std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * y);

            phi_arrs[nbx](i, j, k) =
                radius - std::sqrt(
                             ((x - xc) * (x - xc)) + ((y - yc) * (y - yc)) +
                             ((z - zc) * (z - zc)));
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

void VortexPatch::pre_advance_work()
{
    const auto& time =
        m_sim.time().current_time() + (0.5_rt * m_sim.time().delta_t());

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
                const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
                const amrex::Real xf = problo[0] + (i * dx[0]);
                const amrex::Real yf = problo[1] + (j * dx[1]);
                const amrex::Real zf = problo[2] + (k * dx[2]);
                uf_arrs[nbx](i, j, k) =
                    2.0_rt * std::sin(std::numbers::pi_v<amrex::Real> * xf) *
                    std::sin(std::numbers::pi_v<amrex::Real> * xf) *
                    std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * y) *
                    std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * z) *
                    std::cos(std::numbers::pi_v<amrex::Real> * time / TT);
                vf_arrs[nbx](i, j, k) =
                    -std::sin(std::numbers::pi_v<amrex::Real> * yf) *
                    std::sin(std::numbers::pi_v<amrex::Real> * yf) *
                    std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * x) *
                    std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * z) *
                    std::cos(std::numbers::pi_v<amrex::Real> * time / TT);
                wf_arrs[nbx](i, j, k) =
                    -std::sin(std::numbers::pi_v<amrex::Real> * zf) *
                    std::sin(std::numbers::pi_v<amrex::Real> * zf) *
                    std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * x) *
                    std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * y) *
                    std::cos(std::numbers::pi_v<amrex::Real> * time / TT);
            });
        amrex::Gpu::streamSynchronize();
        u_mac.FillBoundary(geom[lev].periodicity());
        v_mac.FillBoundary(geom[lev].periodicity());
        w_mac.FillBoundary(geom[lev].periodicity());
    }
}

void VortexPatch::post_advance_work()
{
    const auto& time = m_sim.time().current_time();

    const int nlevels = m_sim.repo().num_active_levels();
    const auto& geom = m_sim.mesh().Geom();

    // Overriding the velocity field
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = geom[lev].CellSizeArray();
        const auto& problo = geom[lev].ProbLoArray();
        const amrex::Real TT = m_TT;
        const auto& vel_arrs = m_velocity(lev).arrays();
        amrex::ParallelFor(
            m_velocity(lev), amrex::IntVect(1),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
                const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
                vel_arrs[nbx](i, j, k, 0) =
                    2.0_rt * std::sin(std::numbers::pi_v<amrex::Real> * x) *
                    std::sin(std::numbers::pi_v<amrex::Real> * x) *
                    std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * y) *
                    std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * z) *
                    std::cos(std::numbers::pi_v<amrex::Real> * time / TT);
                vel_arrs[nbx](i, j, k, 1) =
                    -std::sin(std::numbers::pi_v<amrex::Real> * y) *
                    std::sin(std::numbers::pi_v<amrex::Real> * y) *
                    std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * x) *
                    std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * z) *
                    std::cos(std::numbers::pi_v<amrex::Real> * time / TT);
                vel_arrs[nbx](i, j, k, 2) =
                    -std::sin(std::numbers::pi_v<amrex::Real> * z) *
                    std::sin(std::numbers::pi_v<amrex::Real> * z) *
                    std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * x) *
                    std::sin(2.0_rt * std::numbers::pi_v<amrex::Real> * y) *
                    std::cos(std::numbers::pi_v<amrex::Real> * time / TT);
            });
    }
    amrex::Gpu::streamSynchronize();

    m_velocity.fillpatch(time);
}

} // namespace amr_wind
