#include "amr-wind/physics/VortexDipole.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

VortexDipole::VortexDipole(const CFDSim& sim)
    : m_velocity(sim.repo().get_field("velocity"))
    , m_density(sim.repo().get_field("density"))
{
    {
        amrex::ParmParse pp("incflo");
        pp.query("density", m_rho);
    }

    {
        amrex::ParmParse pp("VortexDipole");
        pp.queryarr("left_vortex_location", m_loc_left);
        pp.queryarr("right_vortex_location", m_loc_right);
        pp.query("vortex_core_radius", m_r0);
        pp.queryarr("background_velocity", m_bvel);
    }
}

/** Initialize the velocity and density fields at the beginning of the
 *  simulation.
 */
void VortexDipole::initialize_fields(int level, const amrex::Geometry& geom)
{
    auto& velocity = m_velocity(level);
    auto& density = m_density(level);

    density.setVal(m_rho);

    const auto& problo = geom.ProbLoArray();

    const auto& dx = geom.CellSizeArray();
    const auto& vel_arrs = velocity.arrays();

    const amrex::Real x1 = m_loc_right[0];
    const amrex::Real x2 = m_loc_left[0];

    const amrex::Real z1 = m_loc_right[2];
    const amrex::Real z2 = m_loc_left[2];

    const amrex::Real ub = m_bvel[0];
    const amrex::Real vb = m_bvel[1];
    const amrex::Real wb = m_bvel[2];

    const amrex::Real r0 = m_r0;
    const amrex::Real omegaEmag = m_omegaEmag;

    amrex::ParallelFor(
        velocity, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

            const amrex::Real r1 =
                std::sqrt(((x - x1) * (x - x1)) + ((z - z1) * (z - z1)));

            const amrex::Real r2 =
                std::sqrt(((x - x2) * (x - x2)) + ((z - z2) * (z - z2)));

            vel_arrs[nbx](i, j, k, 0) = ub +
                                        (-0.5_rt * omegaEmag * (z - z1) *
                                         std::exp(-(r1 / r0) * (r1 / r0))) +
                                        (0.5_rt * omegaEmag * (z - z2) *
                                         std::exp(-(r2 / r0) * (r2 / r0)));
            vel_arrs[nbx](i, j, k, 1) = vb;
            vel_arrs[nbx](i, j, k, 2) = wb +
                                        (0.5_rt * omegaEmag * (x - x1) *
                                         std::exp(-(r1 / r0) * (r1 / r0))) +
                                        (-0.5_rt * omegaEmag * (x - x2) *
                                         std::exp(-(r2 / r0) * (r2 / r0)));
        });
    amrex::Gpu::streamSynchronize();
}

} // namespace amr_wind
