#ifndef STOKES_WAVES_K_H
#define STOKES_WAVES_K_H

#include <algorithm>
#include <cmath>
#include <numbers>
#include "AMReX_FArrayBox.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

// Stokes waves theory adapted from
// Fenton, J., Fifth Order Stokes Theory for Steady Waves
// Journal of Waterway, Port, Coastal and Ocean Engineering, 1985, 111, 216-234

// Updated Table 1 from Fenton 1985 found in
// Fenton, J., Nonlinear Wave Theories
// The Sea Vol. 9 Ocean Engineering Science, 1990
// https://johndfenton.com/Papers/Fenton90b-Nonlinear-wave-theories.pdf

// Relevant results are summarized in
// Kinnas S. A., Notes on fifth order gravity wave theory
// https://www.caee.utexas.edu/prof/kinnas/ce358/oenotes/kinnas_stokes11.pdf

namespace amr_wind::ocean_waves::relaxation_zones {

// Compute wavelength as a function of wave period, water depth, and g
AMREX_FORCE_INLINE amrex::Real stokes_wave_length(
    const amrex::Real T,
    const amrex::Real d,
    const amrex::Real H,
    const int order,
    const amrex::Real g,
    const amrex::Real tol,
    const int iter_max)
{
    // Calculate constants and derivatives that do not change with iteration
    const amrex::Real omega = 2.0_rt * std::numbers::pi_v<amrex::Real> / T;
    const amrex::Real depsdk = H / 2.0_rt;

    // Begin Newton-Raphson iterations
    int iter = 0;
    // First guess is first-order
    amrex::Real k = (omega * omega) / g;
    // Cannot skip loop
    amrex::Real f = tol + 1.0_rt;
    while (std::abs(f) > tol && iter < iter_max) {
        // Calculate current constants

        // Exponential definition of S = sech(2kd)
        const amrex::Real S = 2.0_rt * std::exp(2.0_rt * k * d) /
                              (std::exp(4.0_rt * k * d) + 1.0_rt);
        const amrex::Real C = 1.0_rt - S;
        const amrex::Real eps = k * H / 2.0_rt;
        const amrex::Real C0 = std::sqrt(std::tanh(k * d));
        const amrex::Real C2 =
            C0 * (2.0_rt + 7.0_rt * S * S) / (4.0_rt * C * C);
        const amrex::Real numterm_C4 =
            (4.0_rt + (32.0_rt * S) - (116.0_rt * std::pow(S, 2.0_rt)) -
             (400.0_rt * std::pow(S, 3.0_rt)) -
             (71.0_rt * std::pow(S, 4.0_rt)) +
             (146.0_rt * std::pow(S, 5.0_rt)));
        const amrex::Real C4 =
            C0 * numterm_C4 / (32.0_rt * std::pow(C, 5.0_rt));
        // Calculate pure derivates
        const amrex::Real dSdk = -2.0_rt * d * std::sinh(2.0_rt * k * d) /
                                 std::pow(std::cosh(2.0_rt * k * d), 2.0_rt);
        const amrex::Real dCdk = -dSdk;
        const amrex::Real dC0dk =
            d / (2.0_rt * C0 * std::pow(std::cosh(k * d), 2.0_rt));
        // Calculate derivatives with products
        const amrex::Real dC2dk =
            (4.0_rt * std::pow(C, 2.0_rt) *
                 (dC0dk * (2.0_rt + 7.0_rt * std::pow(S, 2.0_rt)) +
                  C0 * 14.0_rt * S * dSdk) -
             C0 * (2.0_rt + 7.0_rt * std::pow(S, 2.0_rt)) * 8.0_rt * C * dCdk) /
            (16.0_rt * std::pow(C, 4.0_rt));
        const amrex::Real dC4dk =
            (32.0_rt * std::pow(C, 5.0_rt) *
                 (dC0dk * numterm_C4 +
                  C0 * (32.0_rt * dSdk - 232.0_rt * S * dSdk -
                        1200.0_rt * std::pow(S, 2.0_rt) * dSdk -
                        284.0_rt * std::pow(S, 3.0_rt) * dSdk +
                        730.0_rt * std::pow(S, 4.0_rt))) -
             C0 * numterm_C4 * 160.0_rt * std::pow(C, 4.0_rt) * dCdk) /
            (1024.0_rt * std::pow(C, 10.0_rt));

        // Calculate derivative for loop convergence
        amrex::Real dfdk =
            (g * std::pow(C0, 2.0_rt)) + (g * k * 2.0_rt * C0 * dC0dk);
        // Add additional terms depending on order
        if (order >= 2) {
            dfdk += (g * (2.0_rt * C0 * std::pow(eps, 2.0_rt) * C2 +
                          std::pow(eps, 4.0_rt) * std::pow(C2, 2.0_rt))) +
                    (g * k *
                     (2.0_rt * dC0dk * std::pow(eps, 2.0_rt) * C2 +
                      2.0_rt * C0 *
                          (2.0_rt * eps * depsdk * C2 +
                           std::pow(eps, 2.0_rt) * dC2dk) +
                      4.0_rt * std::pow(eps, 3.0_rt) * std::pow(C2, 2.0_rt) *
                          depsdk +
                      std::pow(eps, 4.0_rt) * 2.0_rt * C2 * dC2dk));
        }
        if (order >= 4) {
            dfdk +=
                (g * (2.0_rt * std::pow(eps, 4.0_rt) * C0 * C4 +
                      2.0_rt * std::pow(eps, 6.0_rt) * C2 * C4 +
                      std::pow(eps, 8.0_rt) * std::pow(C4, 2.0_rt))) +
                (g * k *
                 (8.0_rt * std::pow(eps, 3.0_rt) * depsdk * C0 * C4 +
                  2.0_rt * std::pow(eps, 4.0_rt) * (C0 * dC4dk + C4 * dC0dk) +
                  12.0_rt * std::pow(eps, 5.0_rt) * depsdk * C2 * C4 +
                  2.0_rt * std::pow(eps, 6.0_rt) * (C2 * dC4dk + C4 * dC2dk) +
                  8.0_rt * std::pow(eps, 7.0_rt) * depsdk *
                      std::pow(C4, 2.0_rt) +
                  std::pow(eps, 8.0_rt) * 2.0_rt * C4 * dC4dk));
        }
        k = k - (f / dfdk);

        iter += 1;
        f = g * k * std::pow(C0, 2.0_rt);
        // Add additional terms depending on order
        if (order >= 2) {
            f += g * k *
                 (2.0_rt * C0 * std::pow(eps, 2.0_rt) * C2 +
                  std::pow(eps, 4.0_rt) * std::pow(C2, 2.0_rt));
        }
        if (order >= 4) {
            f += g * k *
                 (2.0_rt * std::pow(eps, 4.0_rt) * C0 * C4 +
                  2.0_rt * std::pow(eps, 6.0_rt) * C2 * C4 +
                  std::pow(eps, 8.0_rt) * std::pow(C4, 2.0_rt));
        }
        // Subtract omega^2 to measure convergence
        f -= omega * omega;
    }

    if (k < tol) {
        // Return negative wavelength if faulty result
        return -1;
    }
    if (std::isnan(k)) {
        return -2;
    }
    if (iter == iter_max) {
        return -3;
    }
    // Return wavelength calculated from wavenumber
    return 2.0_rt * std::numbers::pi_v<amrex::Real> / k;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void stokes_coefficients(
    int stokes_order,
    amrex::Real wavenumber,
    amrex::Real water_depth,
    amrex::Real& c0,
    amrex::Real& a11,
    amrex::Real& a22,
    amrex::Real& b22,
    amrex::Real& c2,
    amrex::Real& a31,
    amrex::Real& a33,
    amrex::Real& b31,
    amrex::Real& a42,
    amrex::Real& a44,
    amrex::Real& b42,
    amrex::Real& b44,
    amrex::Real& c4,
    amrex::Real& a51,
    amrex::Real& a53,
    amrex::Real& a55,
    amrex::Real& b53,
    amrex::Real& b55)
{

    amrex::Real kd = wavenumber * water_depth;
    kd = amrex::min<amrex::Real>(kd, 50.0_rt * std::numbers::pi_v<amrex::Real>);

    // Exponential definition of S = sech(2kd)
    amrex::Real S =
        2.0_rt * std::exp(2.0_rt * kd) / (std::exp(4.0_rt * kd) + 1.0_rt);
    amrex::Real C = 1.0_rt - S;
    amrex::Real Sh = std::sinh(kd);
    amrex::Real Th = std::tanh(kd);
    // Exponential definition of coth(kd)
    amrex::Real CTh =
        (1.0_rt + std::exp(-2.0_rt * kd)) / (1.0_rt - std::exp(-2 * kd));

    c0 = std::sqrt(Th);
    a11 = 1.0_rt / std::sinh(kd); // Hyperbolic cosecant
    // Second order coefficients
    a22 = 3.0_rt * std::pow(S, 2.0_rt) / (2.0_rt * std::pow(C, 2.0_rt));
    b22 = CTh * (1.0_rt + 2.0_rt * S) / (2.0_rt * C);
    c2 = std::sqrt(Th) * (2.0_rt + 7.0_rt * std::pow(S, 2.0_rt)) /
         (4.0_rt * std::pow(C, 2.0_rt));
    if (stokes_order == 2) {
        return;
    }

    // Third order coefficients
    a31 = (-4.0_rt - 20.0_rt * S + 10.0_rt * std::pow(S, 2.0_rt) -
           13.0_rt * std::pow(S, 3.0_rt)) /
          (8.0_rt * Sh * std::pow(C, 3.0_rt));
    a33 = (-2.0_rt * std::pow(S, 2.0_rt) + 11.0_rt * std::pow(S, 3.0_rt)) /
          (8.0_rt * Sh * std::pow(C, 3.0_rt));
    b31 = -3.0_rt *
          (1.0_rt + 3.0_rt * S + 3.0_rt * std::pow(S, 2.0_rt) +
           2.0_rt * std::pow(S, 3.0_rt)) /
          (8.0_rt * std::pow(C, 3.0_rt));
    if (stokes_order == 3) {
        return;
    }

    // Fourth order coefficients
    a42 = (12.0_rt * S - 14.0_rt * std::pow(S, 2.0_rt) -
           264.0_rt * std::pow(S, 3.0_rt) - 45.0_rt * std::pow(S, 4.0_rt) -
           13.0_rt * std::pow(S, 5.0_rt)) /
          (24.0_rt * std::pow(C, 5.0_rt));
    a44 = (10.0_rt * std::pow(S, 3.0_rt) - 174.0_rt * std::pow(S, 4.0_rt) +
           291.0_rt * std::pow(S, 5.0_rt) + 278.0_rt * std::pow(S, 6.0_rt)) /
          (48.0_rt * (3.0_rt + 2.0_rt * S) * std::pow(C, 5.0_rt));
    b42 = CTh *
          (6.0_rt - 26.0_rt * S - 182.0_rt * std::pow(S, 2.0_rt) -
           204.0_rt * std::pow(S, 3.0_rt) - 25.0_rt * std::pow(S, 4.0_rt) +
           26.0_rt * std::pow(S, 5.0_rt)) /
          (6.0_rt * (3.0_rt + 2.0_rt * S) * std::pow(C, 4.0_rt));
    b44 = CTh *
          (24.0_rt + 92.0_rt * S + 122.0_rt * std::pow(S, 2.0_rt) +
           66.0_rt * std::pow(S, 3.0_rt) + 67.0_rt * std::pow(S, 4.0_rt) +
           34.0_rt * std::pow(S, 5.0_rt)) /
          (24.0_rt * (3.0_rt + 2.0_rt * S) * std::pow(C, 4.0_rt));
    c4 = std::sqrt(Th) *
         (4.0_rt + 32.0_rt * S - 116.0_rt * std::pow(S, 2.0_rt) -
          400.0_rt * std::pow(S, 3.0_rt) - 71.0_rt * std::pow(S, 4.0_rt) +
          146.0_rt * std::pow(S, 5.0_rt)) /
         (32.0_rt * std::pow(C, 5.0_rt));
    if (stokes_order == 4) {
        return;
    }

    // Fifth order coefficients
    a51 = (-1184.0_rt + 32.0_rt * S + 13232.0_rt * std::pow(S, 2.0_rt) +
           21712.0_rt * std::pow(S, 3.0_rt) + 20940.0_rt * std::pow(S, 4.0_rt) +
           12554.0_rt * std::pow(S, 5.0_rt) - 500.0_rt * std::pow(S, 6.0_rt) -
           3341.0_rt * std::pow(S, 7.0_rt) - 670.0_rt * std::pow(S, 8.0_rt)) /
          (64.0_rt * Sh * (3.0_rt + 2.0_rt * S) * (4.0_rt + S) *
           std::pow(C, 6.0_rt));
    a53 = (4.0_rt * S + 105.0_rt * std::pow(S, 2.0_rt) +
           198.0_rt * std::pow(S, 3.0_rt) - 1376.0_rt * std::pow(S, 4.0_rt) -
           1302.0_rt * std::pow(S, 5.0_rt) - 117.0_rt * std::pow(S, 6.0_rt) +
           58.0_rt * std::pow(S, 7.0_rt)) /
          (32.0_rt * Sh * (3.0_rt + 2.0_rt * S) * std::pow(C, 6.0_rt));
    a55 = (-6.0_rt * std::pow(S, 3.0_rt) + 272.0_rt * std::pow(S, 4.0_rt) -
           1552.0_rt * std::pow(S, 5.0_rt) + 852.0_rt * std::pow(S, 6.0_rt) +
           2029.0_rt * std::pow(S, 7.0_rt) + 430.0_rt * std::pow(S, 8.0_rt)) /
          (64.0_rt * Sh * (3.0_rt + 2.0_rt * S) * (4.0_rt + S) *
           std::pow(C, 6.0_rt));
    b53 =
        9.0_rt *
        (132.0_rt + 17.0_rt * S - 2216.0_rt * std::pow(S, 2.0_rt) -
         5897.0_rt * std::pow(S, 3.0_rt) - 6292.0_rt * std::pow(S, 4.0_rt) -
         2687.0_rt * std::pow(S, 5.0_rt) + 194.0_rt * std::pow(S, 6.0_rt) +
         467.0_rt * std::pow(S, 7.0_rt) + 82.0_rt * std::pow(S, 8.0_rt)) /
        (128.0_rt * (3.0_rt + 2.0_rt * S) * (4.0_rt + S) * std::pow(C, 6.0_rt));
    b55 =
        5.0_rt *
        (300.0_rt + 1579.0_rt * S + 3176.0_rt * std::pow(S, 2.0_rt) +
         2949.0_rt * std::pow(S, 3.0_rt) + 1188.0_rt * std::pow(S, 4.0_rt) +
         675.0_rt * std::pow(S, 5.0_rt) + 1326.0_rt * std::pow(S, 6.0_rt) +
         827.0_rt * std::pow(S, 7.0_rt) + 130.0_rt * std::pow(S, 8.0_rt)) /
        (384.0_rt * (3.0_rt + 2.0_rt * S) * (4.0_rt + S) * std::pow(C, 6.0_rt));
    if (stokes_order == 5) {
        return;
    }

    if (stokes_order > 5 || stokes_order < 2) {
        amrex::Abort(
            "invalid stokes order specified. It should be between 2,3,4 or 5 ");
    }
}

// Based on Fenton 1985
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void stokes_waves(
    int stokes_order,
    amrex::Real wavelength,
    amrex::Real water_depth,
    amrex::Real wave_height,
    amrex::Real zsl,
    amrex::Real g,
    amrex::Real x,
    amrex::Real z,
    amrex::Real time,
    amrex::Real phase_offset,
    amrex::Real& eta,
    amrex::Real& u_w,
    amrex::Real& v_w,
    amrex::Real& w_w)
{
    const amrex::Real wavenumber =
        2.0_rt * std::numbers::pi_v<amrex::Real> / wavelength;

    // some parameters
    amrex::Real c0{0.0_rt};
    amrex::Real a11{0.0_rt}, a22{0.0_rt}, b22{0.0_rt}, c2{0.0_rt};
    amrex::Real a31{0.0_rt}, a33{0.0_rt}, b31{0.0_rt}, a42{0.0_rt}, a44{0.0_rt};
    amrex::Real b42{0.0_rt}, b44{0.0_rt}, c4{0.0_rt};
    amrex::Real a51{0.0_rt}, a53{0.0_rt}, a55{0.0_rt}, b53{0.0_rt}, b55{0.0_rt};

    stokes_coefficients(
        stokes_order, wavenumber, water_depth, c0, a11, a22, b22, c2, a31, a33,
        b31, a42, a44, b42, b44, c4, a51, a53, a55, b53, b55);

    const amrex::Real eps = wavenumber * wave_height / 2.0_rt; // Steepness (ka)
    const amrex::Real c =
        (c0 + std::pow(eps, 2.0_rt) * c2 + std::pow(eps, 4.0_rt) * c4) *
        std::sqrt(g / wavenumber);

    const amrex::Real omega = c * wavenumber;
    const amrex::Real phase = (wavenumber * x) - (omega * time) - phase_offset;

    eta = ((eps * std::cos(phase) // first order term
            + std::pow(eps, 2.0_rt) * b22 *
                  std::cos(2.0_rt * phase) // second order term
            + std::pow(eps, 3.0_rt) * b31 *
                  (std::cos(phase) - std::cos(3.0_rt * phase)) +
            std::pow(eps, 4.0_rt) * (b42 * std::cos(2.0_rt * phase) +
                                     b44 * std::cos(4.0_rt * phase)) +
            std::pow(eps, 5.0_rt) * (-(b53 + b55) * std::cos(phase) +
                                     b53 * std::cos(3.0_rt * phase) +
                                     b55 * std::cos(5.0_rt * phase))) /
           wavenumber) +
          zsl;

    // Compute velocities components using Eq.(21) Eq.(23) from Kinnas
    // https://www.caee.utexas.edu/prof/kinnas/ce358/oenotes/kinnas_stokes11.pdf
    // Define coefficients using Eq.(19)

    const int MAX_ORDER = 5;
    amrex::GpuArray<amrex::Real, MAX_ORDER> a;
    if (stokes_order == 2) {
        a[0] = a11;
        a[1] = a22;
    }
    if (stokes_order == 3) {
        a[0] = a11 + ((eps * eps) * a31);
        a[1] = a22;
        a[2] = a33;
    }
    if (stokes_order == 4) {
        a[0] = a11 + ((eps * eps) * a31);
        a[1] = a22 + ((eps * eps) * a42);
        a[2] = a33;
        a[3] = a44;
    }
    if (stokes_order == 5) {
        a[0] = a11 + ((eps * eps) * a31) + (std::pow(eps, 4.0_rt) * a51);
        a[1] = a22 + ((eps * eps) * a42);
        a[2] = a33 + ((eps * eps) * a53);
        a[3] = a44;
        a[4] = a55;
    }

    u_w = 0.0_rt;
    v_w = 0.0_rt;
    w_w = 0.0_rt;
    for (int n = 1; n <= stokes_order; ++n) {
        // Upper bound for deep water case
        // This ensure finite values of velocity for large kd's
        amrex::Real nkdz = n * wavenumber * (water_depth + (z - zsl));
        nkdz = amrex::min<amrex::Real>(
            nkdz, 50.0_rt * std::numbers::pi_v<amrex::Real>);
        u_w += std::pow(eps, static_cast<amrex::Real>(n)) *
               static_cast<amrex::Real>(n) * a[n - 1] * std::cosh(nkdz) *
               std::cos(static_cast<amrex::Real>(n) * phase);
        w_w += std::pow(eps, static_cast<amrex::Real>(n)) *
               static_cast<amrex::Real>(n) * a[n - 1] * std::sinh(nkdz) *
               std::sin(static_cast<amrex::Real>(n) * phase);
    }
    u_w *= (c0 * std::sqrt(g / wavenumber));
    w_w *= (c0 * std::sqrt(g / wavenumber));
}

} // namespace amr_wind::ocean_waves::relaxation_zones
#endif
