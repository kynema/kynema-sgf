#ifndef MOSD_H
#define MOSD_H

#include <numbers>
#include "AMReX_AmrCore.H"
#include "AMReX.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;
namespace amr_wind {
struct MOSD
{
    /*
     * A dynamic wall model that calculates the stress from wave to wind
     * based on geometric information of the wave.
     */

    amrex::Real amplitude{0.05_rt};
    amrex::Real wavenumber{4.0_rt};
    amrex::Real omega{0.8_rt};
    amrex::Real time;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real get_dyn_tau(
        const amrex::Real u_dx,
        const amrex::Real v_dx,
        const amrex::Real xc,
        const amrex::Real unit_nor) const
    {

        // Building the wave surface, gradients, wave velocities and unit normal
        const amrex::Real dx_eta_wave =
            -amplitude * wavenumber *
            std::sin((wavenumber * xc) - (omega * time));
        const amrex::Real dy_eta_wave = 0;
        const amrex::Real dt_eta_wave =
            amplitude * omega * std::sin((wavenumber * xc) - (omega * time));
        const amrex::Real grad_eta_wave = std::sqrt(
            (dx_eta_wave * dx_eta_wave) + (dy_eta_wave * dy_eta_wave));
        const amrex::Real Cx_wave = -dt_eta_wave * dx_eta_wave *
                                    (1.0_rt / (grad_eta_wave * grad_eta_wave));
        const amrex::Real Cy_wave = -dt_eta_wave * dy_eta_wave *
                                    (1.0_rt / (grad_eta_wave * grad_eta_wave));
        const amrex::Real n_x = dx_eta_wave / grad_eta_wave;
        const amrex::Real n_y = dy_eta_wave / grad_eta_wave;

        // Calculating the relative velocity, heaviside function
        const amrex::Real u_r = u_dx - Cx_wave;
        const amrex::Real v_r = v_dx - Cy_wave;
        const amrex::Real ur_mag =
            std::sqrt((u_r * u_r * n_x * n_x) + (v_r * v_r * n_y * n_y));
        const amrex::Real Heavi_arg =
            ((u_r * dx_eta_wave) + (v_r * dy_eta_wave));
        const amrex::Real Heavi =
            (Heavi_arg + std::abs(Heavi_arg)) / (2.0_rt * Heavi_arg);

        // Calculating the magnitude of the stress
        return (1.0_rt / std::numbers::pi_v<amrex::Real>)*ur_mag * ur_mag *
               grad_eta_wave * grad_eta_wave * Heavi *
               (unit_nor == 0 ? n_x : n_y);
    }
};

} // namespace amr_wind
#endif /* MOSD_H */
