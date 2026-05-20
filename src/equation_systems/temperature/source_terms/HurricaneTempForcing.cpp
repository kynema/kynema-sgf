#include "src/equation_systems/temperature/source_terms/HurricaneTempForcing.H"
#include "src/CFDSim.H"
#include "src/utilities/trig_ops.H"
#include "src/core/vs/vstraits.H"
#include "src/wind_energy/ABL.H"
#include "src/core/FieldUtils.H"
#include "src/utilities/linear_interpolation.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::temperature {

HurricaneTempForcing::HurricaneTempForcing(const CFDSim& sim)
    : m_mesh(sim.mesh())
{

    const auto& abl = sim.physics_manager().get<kynema_sgf::ABL>();
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

void HurricaneTempForcing::operator()(
    const int lev, const FieldState /*fstate*/, amrex::MultiFab& src_term) const
{
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    const amrex::Real dTdR = m_dTdR;
    const amrex::Real dTzh = m_dTzh;

    const int idir = m_axis;
    const amrex::Real* heights = m_vel_ht.data();
    const amrex::Real* heights_end = m_vel_ht.end();
    const amrex::Real* vals = m_vel_vals.data();

    auto const& src_arrs = src_term.arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), 1,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
            amrex::IntVect iv(i, j, k);
            const amrex::Real ht =
                problo[idir] + ((iv[idir] + 0.5_rt) * dx[idir]);

            const amrex::Real dTdR_z = dTdR * (dTzh - ht) / dTzh;
            const amrex::Real vmean = kynema_sgf::interp::linear(
                heights, heights_end, vals, ht, 3, 1);

            src_arrs[nbx](i, j, k) -= vmean * dTdR_z;
        });
}

} // namespace kynema_sgf::pde::temperature
