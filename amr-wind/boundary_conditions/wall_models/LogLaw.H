#ifndef LOGLAW_H
#define LOGLAW_H

#include "AMReX_AmrCore.H"
#include "AMReX.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;
namespace amr_wind {
struct LogLaw
{
    /*
     * A simple wall model that sets the wall-shear stress
     * based on computing u_tau given the horizontal velocity
     * magnitude at a zref. This is akin to an explicit non-linear
     * Robin boundary condition at the wall.
     */

    // Log law constants from Lee & Moser 2015
    // https://doi.org/10.1017/jfm.2015.268.
    amrex::Real B{4.27_rt};
    amrex::Real kappa{0.384_rt};
    int max_iters = 25; // Max iterations for u_tau Newton-Raphson solve
    // Reference height for log law
    amrex::Real zref;
    int ref_index{0};
    amrex::Real nu; // molecular viscosity
    // u_tau state variable, gets updated in update_utau depending on
    // the type of wall model used
    amrex::Real utau_mean{1.0_rt};
    amrex::Real wspd_mean; // mean horizontal velocity magnitude

    void update_utau_mean() { utau_mean = get_utau(wspd_mean); }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    get_utau(amrex::Real wspd) const
    {
        amrex::Real utau_iter = -1.0_rt;
        amrex::Real wspd_pred;
        amrex::Real wspd_deriv;
        amrex::Real zplus;
        amrex::Real utau = utau_mean;
        int iter = 0;
        while ((std::abs(utau_iter - utau) > 1.0e-5_rt) && iter <= max_iters) {
            utau_iter = utau;
            zplus = zref * utau / nu;
            // Get wspd for a given utau from log-law
            wspd_pred = utau * (std::log(zplus) / kappa + B);
            wspd_deriv =
                ((1.0_rt + std::log(zplus)) / kappa) + B; // d(wspd)/d(utau)
            utau = utau -
                   ((wspd_pred - wspd) / wspd_deriv); // Newton-Raphson update
            ++iter;
        }
        if (iter == max_iters) {
            amrex::Abort();
        }
        return utau;
    }
};
} /* namespace amr_wind */

#endif /* LogLaw_H */
