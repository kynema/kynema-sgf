#include "amr-wind/equation_systems/icns/source_terms/RayleighDamping.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::icns {

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
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{

    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& probhi = m_mesh.Geom(lev).ProbHiArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    const amrex::Real tau = m_tau;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> ref_vel{
        m_ref_vel[0], m_ref_vel[1], m_ref_vel[2]};
    const auto& vel =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_array(mfi);

    // Constants used to determine the fringe region coefficient
    const amrex::Real dRD = m_dRD;
    const amrex::Real dFull = m_dFull;

    // Which coordinate directions to force
    const amrex::Real fx = m_fcoord[0];
    const amrex::Real fy = m_fcoord[1];
    const amrex::Real fz = m_fcoord[2];

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
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
        src_term(i, j, k, 0) +=
            fx * coeff * (ref_vel[0] - vel(i, j, k, 0)) / tau;
        src_term(i, j, k, 1) +=
            fy * coeff * (ref_vel[1] - vel(i, j, k, 1)) / tau;
        src_term(i, j, k, 2) +=
            fz * coeff * (ref_vel[2] - vel(i, j, k, 2)) / tau;
    });
}

} // namespace amr_wind::pde::icns
