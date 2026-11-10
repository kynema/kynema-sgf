#include "src/physics/ForestCanopy.H"
#include "AMReX_ParmParse.H"
#include <cmath>

namespace kynema_sgf::forestdrag {

CanopyTurbulence parse_canopy_turbulence()
{
    CanopyTurbulence canopy;
    amrex::ParmParse pp("ForestDrag");
    pp.query("canopy_tke", canopy.m_enabled);
    pp.query("canopy_beta_p", canopy.m_beta_p);
    pp.query("canopy_beta_d", canopy.m_beta_d);
    pp.query("canopy_length_alpha", canopy.m_length_alpha);
    if (canopy.m_beta_p < 0.0_rt || canopy.m_beta_d < 0.0_rt ||
        canopy.m_length_alpha < 0.0_rt) {
        amrex::Abort(
            "ForestDrag: canopy_beta_p, canopy_beta_d and canopy_length_alpha "
            "must not be negative");
    }
    return canopy;
}

void add_canopy_tke_source(
    amrex::MultiFab& src_term,
    const amrex::MultiFab& forest_drag,
    const amrex::MultiFab& velocity,
    const amrex::MultiFab& tke,
    const amrex::Real beta_p,
    const amrex::Real beta_d,
    const amrex::Real dt)
{
    auto const& src_arrs = src_term.arrays();
    auto const& forest_arrs = forest_drag.const_arrays();
    auto const& vel_arrs = velocity.const_arrays();
    auto const& tke_arrs = tke.const_arrays();
    amrex::ParallelFor(
        src_term, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real fd = forest_arrs[nbx](i, j, k);
            if (fd > 0.0_rt) {
                const amrex::Real ux = vel_arrs[nbx](i, j, k, 0);
                const amrex::Real uy = vel_arrs[nbx](i, j, k, 1);
                const amrex::Real uz = vel_arrs[nbx](i, j, k, 2);
                const amrex::Real ws =
                    std::sqrt((ux * ux) + (uy * uy) + (uz * uz));
                // expm1 is in (-1, 0]: the sink never exceeds k / dt
                src_arrs[nbx](i, j, k) +=
                    (beta_p * fd * ws * ws * ws) +
                    (tke_arrs[nbx](i, j, k) *
                     std::expm1(-beta_d * fd * ws * dt) / dt);
            }
        });
}

void limit_canopy_length_scale(
    const amrex::MultiFab& forest_drag,
    const amrex::Real alpha,
    amrex::MultiFab& turb_lscale,
    amrex::MultiFab& mu_turb,
    amrex::MultiFab& shear_prod,
    amrex::MultiFab& buoy_prod)
{
    auto const& forest_arrs = forest_drag.const_arrays();
    auto const& tlscale_arrs = turb_lscale.arrays();
    auto const& mu_arrs = mu_turb.arrays();
    auto const& shear_prod_arrs = shear_prod.arrays();
    auto const& buoy_prod_arrs = buoy_prod.arrays();
    amrex::ParallelFor(
        mu_turb, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real fd = forest_arrs[nbx](i, j, k);
            const amrex::Real lscale = tlscale_arrs[nbx](i, j, k);
            // alpha / fd < lscale, written without dividing by a small fd
            const amrex::Real lscale_fd = lscale * fd;
            if (fd > 0.0_rt && alpha < lscale_fd) {
                const amrex::Real ratio = alpha / lscale_fd;
                tlscale_arrs[nbx](i, j, k) = alpha / fd;
                mu_arrs[nbx](i, j, k) *= ratio;
                shear_prod_arrs[nbx](i, j, k) *= ratio;
                buoy_prod_arrs[nbx](i, j, k) *= ratio;
            }
        });
}

} // namespace kynema_sgf::forestdrag
