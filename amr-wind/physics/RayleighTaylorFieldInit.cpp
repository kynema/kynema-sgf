#include <cmath>
#include <numbers>

#include "amr-wind/physics/RayleighTaylorFieldInit.H"
#include "AMReX_Gpu.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

RayleighTaylorFieldInit::RayleighTaylorFieldInit()
{
    amrex::ParmParse pp("RayleighTaylor");
    pp.query("rho_lo", m_rho_lo);
    pp.query("rho_hi", m_rho_hi);
}

void RayleighTaylorFieldInit::operator()(
    const amrex::Box& vbx,
    const amrex::Geometry& geom,
    const amrex::Array4<amrex::Real>& density) const
{
    const auto& dx = geom.CellSizeArray();
    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();

    const amrex::Real rho_1 = m_rho_lo;
    const amrex::Real rho_2 = m_rho_hi;
    const amrex::Real splitx = 0.5_rt * (problo[0] + probhi[0]);
    const amrex::Real splity = 0.5_rt * (problo[1] + probhi[1]);
    const amrex::Real L_x = probhi[0] - problo[0];

    amrex::ParallelFor(vbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
        const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
        const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
        const amrex::Real r2d = amrex::min<amrex::Real>(
            std::hypot((x - splitx), (y - splity)), 0.5_rt * L_x);
        const amrex::Real pertheight =
            0.5_rt -
            (0.01_rt *
             std::cos(2.0_rt * std::numbers::pi_v<amrex::Real> * r2d / L_x));

        density(i, j, k) =
            rho_1 + (((rho_2 - rho_1) / 2.0_rt) *
                     (1.0_rt + std::tanh((z - pertheight) / 0.005_rt)));
    });
}

} // namespace amr_wind
