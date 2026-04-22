#ifndef INTEGRALS_H
#define INTEGRALS_H
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::utils {

/** Trapezoid integration
 *
 *  Prefer exact evaluation of the integral over using this function
 */
template <typename Function>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real trapz(
    const amrex::Real xa, const amrex::Real xb, const int n, const Function& f)
{

    const amrex::Real dx = (xb - xa) / n;

    // first and last terms
    amrex::Real sum = 0.5_rt * (f(xa) + f(xb));

    for (int i = 1; i < n; i++) {
        sum += f(xa + (i * dx));
    }
    return dx * sum;
}

} // namespace amr_wind::utils
#endif
