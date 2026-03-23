#include "amr-wind/equation_systems/temperature/source_terms/HurricaneTempForcing.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/core/vs/vstraits.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/utilities/linear_interpolation.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::temperature {

HurricaneTempForcing::HurricaneTempForcing(const CFDSim& sim)
    : m_mesh(sim.mesh())
{

    const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
    // NO need to re-register the hurricane forcing
    abl.register_hurricane_temp_forcing(this);
    // Read the Hurricane Temperature Forcing
    {
        amrex::ParmParse pp("HurricaneTempForcing");
        pp.query("radial_decay", m_dTdR);
        pp.query("radial_decay_zero_height", m_dTzh);

        mean_velocity_init(abl.abl_statistics().vel_profile_coarse());
    }
}

HurricaneTempForcing::~HurricaneTempForcing() = default;

void HurricaneTempForcing::operator()(
    const int lev,
    const amrex::MFIter& /*mfi*/,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const

{
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    const amrex::Real dTdR = m_dTdR;
    const amrex::Real dTzh = m_dTzh;

    // Mean velocity profile used to compute background hurricane forcing term

    const int idir = m_axis;
    const amrex::Real* heights = m_vel_ht.data();
    const amrex::Real* heights_end = m_vel_ht.end();
    const amrex::Real* vals = m_vel_vals.data();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        amrex::IntVect iv(i, j, k);
        const amrex::Real ht = problo[idir] + ((iv[idir] + 0.5_rt) * dx[idir]);

        /*const amrex::Real umean =
            vals[3 * il + 0] + ((vals[3 * ir + 0] - vals[3 * il + 0]) /
                                (heights[ir] - heights[il])) *
                                   (ht - heights[il]);
        */
        const amrex::Real dTdR_z = dTdR * (dTzh - ht) / dTzh;
        const amrex::Real vmean =
            amr_wind::interp::linear(heights, heights_end, vals, ht, 3, 1);

        src_term(i, j, k) -= vmean * dTdR_z;
    });
}

void HurricaneTempForcing::mean_velocity_init(const VelPlaneAveraging& vavg)
{
    m_axis = vavg.axis();

    // The implementation depends the assumption that the ABL statistics
    // class computes statistics at the cell-centeres only on level 0. If
    // this assumption changes in future, the implementation will break...
    // so put in a check here to catch this.
    AMREX_ALWAYS_ASSERT(
        m_mesh.Geom(0).Domain().length(m_axis) ==
        static_cast<int>(vavg.line_centroids().size()));

    m_vel_ht.resize(vavg.line_centroids().size());
    m_vel_vals.resize(vavg.line_average().size());

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, vavg.line_centroids().begin(),
        vavg.line_centroids().end(), m_vel_ht.begin());

    mean_velocity_update(vavg);
}

void HurricaneTempForcing::mean_velocity_update(const VelPlaneAveraging& vavg)
{
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, vavg.line_average().begin(),
        vavg.line_average().end(), m_vel_vals.begin());
}

} // namespace amr_wind::pde::temperature
