#include <limits>
#include "amr-wind/wind_energy/MOData.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

/*
 * van der Laan, P., Kelly, M. C., & Sørensen, N. N. (2017). A new k-epsilon
 * model consistent with Monin-Obukhov similarity theory. Wind Energy,
 * 20(3), 479–489. https://doi.org/10.1002/we.2017
 *
 * Consistent with Dyer (1974) formulation from page 57, Chapter 2, Modeling
 * the vertical ABL structure in Modelling of Atmospheric Flow Fields,
 * Demetri P Lalas and Corrado F Ratto, January 1996,
 * https://doi.org/10.1142/2975.
 */

amrex::Real MOData::calc_psi_m(
    const amrex::Real zeta, const amrex::Real beta_m, const amrex::Real gamma_m)
{
    if (zeta > 0) {
        return -gamma_m * zeta;
    }
    const amrex::Real x = std::sqrt(std::sqrt(1.0_rt - (beta_m * zeta)));
    return (2.0_rt * std::log(0.5_rt * (1.0_rt + x))) +
           std::log(0.5_rt * (1.0_rt + x * x)) - (2.0_rt * std::atan(x)) +
           utils::half_pi();
}

amrex::Real MOData::calc_psi_m(amrex::Real zeta) const
{
    return calc_psi_m(zeta, beta_m, gamma_m);
}

amrex::Real MOData::calc_psi_h(
    const amrex::Real zeta, const amrex::Real beta_h, const amrex::Real gamma_h)
{
    if (zeta > 0) {
        return -gamma_h * zeta;
    }
    const amrex::Real x = std::sqrt(1 - (beta_h * zeta));
    return 2.0_rt * std::log(0.5_rt * (1 + x));
}

amrex::Real MOData::calc_psi_h(amrex::Real zeta) const
{
    return calc_psi_h(zeta, beta_h, gamma_h);
}

void MOData::update_fluxes(int max_iters)
{
    constexpr amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();
    amrex::Real zeta = 0.0_rt;
    amrex::Real utau_iter = 0.0_rt;

    // Initialize variables
    amrex::Real psi_m = 0.0_rt;
    amrex::Real psi_h = 0.0_rt;
    utau = kappa * vmag_mean / (std::log(zref / z0));

    int iter = 0;
    do {
        utau_iter = utau;
        switch (alg_type) {
        case ThetaCalcType::HEAT_FLUX:
            surf_temp = (surf_temp_flux * (std::log(zref / z0t) - psi_h) /
                         (utau * kappa)) +
                        theta_mean;
            break;

        case ThetaCalcType::SURFACE_TEMPERATURE:
            surf_temp_flux = -(theta_mean - surf_temp) * utau * kappa /
                             (std::log(zref / z0t) - psi_h);
            break;
        }

        if (std::abs(surf_temp_flux) > eps) {
            // Stable and unstable ABL conditions
            obukhov_len = -utau * utau * utau * theta_mean /
                          (kappa * gravity * surf_temp_flux);
            zeta = zref / obukhov_len;
        } else {
            // Neutral conditions
            obukhov_len = std::numeric_limits<amrex::Real>::max();
            zeta = 0.0_rt;
        }
        psi_m = calc_psi_m(zeta);
        psi_h = calc_psi_h(zeta);
        utau = kappa * vmag_mean / (std::log(zref / z0) - psi_m);
        ++iter;
    } while ((std::abs(utau_iter - utau) > 1.0e-5_rt) && iter <= max_iters);

    if (iter >= max_iters) {
        amrex::Print()
            << "MOData::update_fluxes: Convergence criteria not met after "
            << max_iters << " iterations\nObuhov length = " << obukhov_len
            << " zeta = " << zeta << "\npsi_m = " << psi_m
            << " psi_h = " << psi_h << "\nutau = " << utau
            << " Tsurf = " << surf_temp << " q = " << surf_temp_flux << '\n';
    }
}

} // namespace amr_wind
