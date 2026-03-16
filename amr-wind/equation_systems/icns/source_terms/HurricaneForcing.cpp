#include "amr-wind/equation_systems/icns/source_terms/HurricaneForcing.H"
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

namespace amr_wind::pde::icns {

HurricaneForcing::HurricaneForcing(const CFDSim& sim) : m_mesh(sim.mesh())
{

    const auto& abl = sim.physics_manager().get<amr_wind::ABL>();
    abl.register_hurricane_forcing(this);

    {
        // Read the rotational time period (in seconds)
        amrex::ParmParse pp("CoriolisForcing");
        amrex::Real rot_time_period = 86400.0_rt;
        pp.query("rotational_time_period", rot_time_period);
        amrex::Real latitude = 90.0_rt;
        pp.query("latitude", latitude);
        m_coriolis_factor = (2.0_rt * utils::two_pi() / rot_time_period) *
                            std::sin(utils::radians(latitude));
        amrex::Print() << "Geostrophic forcing: Coriolis factor = "
                       << m_coriolis_factor << '\n';
    }

    {
        // Read the geostrophic wind speed vector (in m/s)
        amrex::ParmParse pp("HurricaneForcing");
        pp.query("gradient_wind", m_V);
        pp.query("gradient_wind_radial_decay", m_dVdR);
        pp.query("gradient_wind_zero_height", m_Vzh);
        pp.query("eyewall_radial_distance", m_R);
    }

    mean_velocity_init(abl.abl_statistics().vel_profile_coarse());
}

HurricaneForcing::~HurricaneForcing() = default;

void HurricaneForcing::operator()(
    const int lev,
    const amrex::MFIter& /*mfi*/,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{

    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();

    const amrex::Real R = m_R;
    const amrex::Real V = m_V;
    const amrex::Real dVdR = m_dVdR;
    const amrex::Real Vzh = m_Vzh;
    const amrex::Real f = m_coriolis_factor;

    // Mean velocity profile used to compute background hurricane forcing term

    const int idir = m_axis;
    const amrex::Real* heights = m_vel_ht.data();
    const amrex::Real* heights_end = m_vel_ht.end();
    const amrex::Real* vals = m_vel_vals.data();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        amrex::IntVect iv(i, j, k);
        const amrex::Real ht = problo[idir] + ((iv[idir] + 0.5_rt) * dx[idir]);

        const amrex::Real umean =
            amr_wind::interp::linear(heights, heights_end, vals, ht, 3, 0);
        const amrex::Real vmean =
            amr_wind::interp::linear(heights, heights_end, vals, ht, 3, 1);

        // Gradient velocities varying with height
        const amrex::Real V_z = V * (Vzh - ht) / Vzh;
        const amrex::Real dVdR_z = dVdR * (Vzh - ht) / Vzh;
        // Compute the LES terms as presented in George Bryan's paper
        const amrex::Real M1LES =
            (umean * umean / R) + (vmean * V_z / R) - (f * V_z + V_z * V_z / R);
        const amrex::Real M2LES = (-umean * dVdR_z) - (umean * V_z / R);

        src_term(i, j, k, 0) += M1LES;
        src_term(i, j, k, 1) += M2LES;
        src_term(i, j, k, 2) += 0.0_rt;
    });
}

void HurricaneForcing::mean_velocity_init(const VelPlaneAveraging& vavg)
{
    m_axis = vavg.axis();

    // The implementation depends the assumption that the ABL statistics class
    // computes statistics at the cell-centeres only on level 0. If this
    // assumption changes in future, the implementation will break... so put in
    // a check here to catch this.
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

void HurricaneForcing::mean_velocity_update(const VelPlaneAveraging& vavg)
{
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, vavg.line_average().begin(),
        vavg.line_average().end(), m_vel_vals.begin());
}

} // namespace amr_wind::pde::icns
