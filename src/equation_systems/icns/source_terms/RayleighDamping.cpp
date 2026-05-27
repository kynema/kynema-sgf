#include "src/equation_systems/icns/source_terms/RayleighDamping.H"
#include "src/CFDSim.H"
#include "src/utilities/trig_ops.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::icns {

RayleighDamping::RayleighDamping(const CFDSim& sim)
    : m_mesh(sim.mesh()), m_velocity(sim.repo().get_field("velocity"))
{
    // Read the Rayleigh Damping Layer parameters
    amrex::ParmParse pp("RayleighDamping");
    pp.get("time_scale", m_tau);
    // Length where damping coefficient depends on spatial position
    // In sloped region, coefficient goes from 1 to 0
    pp.get("length_sloped_damping", m_dRD);
    // Length where damping coefficient is set to 1
    pp.get("length_complete_damping", m_dFull);
    // Total damping length is m_dRD + m_dFull. Total length is not read in.
    pp.getarr("reference_velocity", m_ref_vel);

    // Which coordinate directions to force
    pp.queryarr("force_coord_directions", m_fcoord);

    // Based upon Allaerts & Meyers (JFM, 2017) and Durran & Klemp (AMS, 1983)
}

RayleighDamping::~RayleighDamping() = default;

void RayleighDamping::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& probhi = m_mesh.Geom(lev).ProbHiArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    const amrex::Real tau = m_tau;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> ref_vel{
        m_ref_vel[0], m_ref_vel[1], m_ref_vel[2]};

    const amrex::Real dRD = m_dRD;
    const amrex::Real dFull = m_dFull;

    const amrex::Real fx = m_fcoord[0];
    const amrex::Real fy = m_fcoord[1];
    const amrex::Real fz = m_fcoord[2];

    auto const& src_arrs = src_term.arrays();
    auto const& vel_arrs =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
            amrex::Real coeff = 0.0_rt;
            const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

            if (probhi[2] - z > dRD + dFull) {
                coeff = 0.0_rt;
            } else if (probhi[2] - z > dFull) {
                coeff = (0.5_rt * std::cos(
                                      std::numbers::pi_v<amrex::Real> *
                                      (probhi[2] - dFull - z) / dRD)) +
                        0.5_rt;
            } else {
                coeff = 1.0_rt;
            }
            src_arrs[nbx](i, j, k, n) +=
                ((n == 0) ? fx : ((n == 1) ? fy : fz)) * coeff *
                (ref_vel[n] - vel_arrs[nbx](i, j, k, n)) / tau;
        });
}

} // namespace kynema_sgf::pde::icns
