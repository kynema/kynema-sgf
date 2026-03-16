#ifndef ANALYTICAL_FUNCTION_H
#define ANALYTICAL_FUNCTION_H

#include "aw_test_utils/AmrexTest.H"
#include "aw_test_utils/MeshTest.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests::analytical_function {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real phi_eval(
    int degree,
    const amrex::Real* coeff,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    amrex::Real phi = 0.0_rt;

    int ind = 0;
    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j) {
            for (int k = 0; k <= degree; ++k) {
                phi += coeff[ind++] * static_cast<amrex::Real>(std::pow(x, i)) *
                       static_cast<amrex::Real>(std::pow(y, j)) *
                       static_cast<amrex::Real>(std::pow(z, k));
            }
        }
    }

    return phi;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real filter_eval(
    int degree,
    const amrex::Real* coeff,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z,
    amrex::Real dx,
    amrex::Real dy,
    amrex::Real dz)
{
    amrex::Real filt_phi = 0.0_rt;
    amrex::Real filt_phi_x = 0.0_rt;
    amrex::Real filt_phi_y = 0.0_rt;
    amrex::Real filt_phi_z = 0.0_rt;
    int ind = 0;

    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j) {
            for (int k = 0; k <= degree; ++k) {
                amrex::Real C = coeff[ind++];
                filt_phi_x +=
                    1.0_rt / (2.0_rt * dx) * C *
                    (static_cast<amrex::Real>(std::pow(x + dx, i + 1)) -
                     static_cast<amrex::Real>(std::pow(x - dx, i + 1))) /
                    (i + 1) * static_cast<amrex::Real>(std::pow(y, j)) *
                    static_cast<amrex::Real>(std::pow(z, k));
                filt_phi_y +=
                    1.0_rt / (2.0_rt * dy) * C *
                    static_cast<amrex::Real>(std::pow(x, i)) *
                    (static_cast<amrex::Real>(std::pow(y + dy, j + 1)) -
                     static_cast<amrex::Real>(std::pow(y - dy, j + 1))) /
                    (j + 1) * static_cast<amrex::Real>(std::pow(z, k));
                filt_phi_z +=
                    1.0_rt / (2.0_rt * dz) * C *
                    static_cast<amrex::Real>(std::pow(x, i)) *
                    static_cast<amrex::Real>(std::pow(y, j)) *
                    (static_cast<amrex::Real>(std::pow(z + dz, k + 1)) -
                     static_cast<amrex::Real>(std::pow(z - dz, k + 1))) /
                    (k + 1);
            }
        }
    }
    filt_phi = 1.0_rt / 3.0_rt * (filt_phi_x + filt_phi_y + filt_phi_z);
    return filt_phi;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real dphidx_eval(
    int degree,
    const amrex::Real* coeff,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    amrex::Real dphidx = 0.0_rt;

    int ind = 0;
    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j) {
            for (int k = 0; k <= degree; ++k) {
                dphidx += i * coeff[ind++] *
                          static_cast<amrex::Real>(
                              std::pow(x, amrex::max(i - 1, 0))) *
                          static_cast<amrex::Real>(std::pow(y, j)) *
                          static_cast<amrex::Real>(std::pow(z, k));
            }
        }
    }

    return dphidx;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real dphidy_eval(
    int degree,
    const amrex::Real* coeff,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    amrex::Real dphidy = 0.0_rt;

    int ind = 0;
    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j) {
            for (int k = 0; k <= degree; ++k) {
                dphidy += j * coeff[ind++] *
                          static_cast<amrex::Real>(std::pow(x, i)) *
                          static_cast<amrex::Real>(
                              std::pow(y, amrex::max(j - 1, 0))) *
                          static_cast<amrex::Real>(std::pow(z, k));
            }
        }
    }

    return dphidy;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real dphidz_eval(
    int degree,
    const amrex::Real* coeff,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    amrex::Real dphidz = 0.0_rt;

    int ind = 0;
    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j) {
            for (int k = 0; k <= degree; ++k) {
                dphidz +=
                    k * coeff[ind++] *
                    static_cast<amrex::Real>(std::pow(x, i)) *
                    static_cast<amrex::Real>(std::pow(y, j)) *
                    static_cast<amrex::Real>(std::pow(z, amrex::max(k - 1, 0)));
            }
        }
    }

    return dphidz;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real laplacian(
    int degree,
    const amrex::Real* coeffx,
    const amrex::Real* coeffy,
    const amrex::Real* coeffz,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    amrex::Real lap = 0.0_rt;

    int ind = 0;
    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j) {
            for (int k = 0; k <= degree; ++k) {
                lap += i * amrex::max(i - 1, 0) * coeffx[ind] *
                       static_cast<amrex::Real>(
                           std::pow(x, amrex::max(i - 2, 0))) *
                       static_cast<amrex::Real>(std::pow(y, j)) *
                       static_cast<amrex::Real>(std::pow(z, k));
                lap += j * amrex::max(j - 1, 0) * coeffy[ind] *
                       static_cast<amrex::Real>(std::pow(x, i)) *
                       static_cast<amrex::Real>(
                           std::pow(y, amrex::max(j - 2, 0))) *
                       static_cast<amrex::Real>(std::pow(z, k));
                lap +=
                    k * amrex::max(k - 1, 0) * coeffz[ind] *
                    static_cast<amrex::Real>(std::pow(x, i)) *
                    static_cast<amrex::Real>(std::pow(y, j)) *
                    static_cast<amrex::Real>(std::pow(z, amrex::max(k - 2, 0)));

                ++ind;
            }
        }
    }

    return lap;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real divergence(
    int degree,
    const amrex::Real* coeffx,
    const amrex::Real* coeffy,
    const amrex::Real* coeffz,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    amrex::Real div = 0.0_rt;

    int ind = 0;
    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j) {
            for (int k = 0; k <= degree; ++k) {
                div += i * coeffx[ind] *
                       static_cast<amrex::Real>(
                           std::pow(x, amrex::max(i - 1, 0))) *
                       static_cast<amrex::Real>(std::pow(y, j)) *
                       static_cast<amrex::Real>(std::pow(z, k));
                div += j * coeffy[ind] *
                       static_cast<amrex::Real>(std::pow(x, i)) *
                       static_cast<amrex::Real>(
                           std::pow(y, amrex::max(j - 1, 0))) *
                       static_cast<amrex::Real>(std::pow(z, k));
                div +=
                    k * coeffz[ind] * static_cast<amrex::Real>(std::pow(x, i)) *
                    static_cast<amrex::Real>(std::pow(y, j)) *
                    static_cast<amrex::Real>(std::pow(z, amrex::max(k - 1, 0)));

                ++ind;
            }
        }
    }

    return div;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real strainrate(
    int degree,
    const amrex::Real* coeffx,
    const amrex::Real* coeffy,
    const amrex::Real* coeffz,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    const amrex::Real ux = dphidx_eval(degree, coeffx, x, y, z);
    const amrex::Real uy = dphidy_eval(degree, coeffx, x, y, z);
    const amrex::Real uz = dphidz_eval(degree, coeffx, x, y, z);

    const amrex::Real vx = dphidx_eval(degree, coeffy, x, y, z);
    const amrex::Real vy = dphidy_eval(degree, coeffy, x, y, z);
    const amrex::Real vz = dphidz_eval(degree, coeffy, x, y, z);

    const amrex::Real wx = dphidx_eval(degree, coeffz, x, y, z);
    const amrex::Real wy = dphidy_eval(degree, coeffz, x, y, z);
    const amrex::Real wz = dphidz_eval(degree, coeffz, x, y, z);

    return std::sqrt(
        (2.0_rt * ux * ux) + (2.0_rt * vy * vy) + (2.0_rt * wz * wz) +
        ((uy + vx) * (uy + vx)) + ((vz + wy) * (vz + wy)) +
        ((wx + uz) * (wx + uz)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real nonlinearsum_11(
    int degree,
    const amrex::Real* coeffx,
    const amrex::Real* coeffy,
    const amrex::Real* coeffz,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    const amrex::Real ux = dphidx_eval(degree, coeffx, x, y, z);
    const amrex::Real uy = dphidy_eval(degree, coeffx, x, y, z);
    const amrex::Real uz = dphidz_eval(degree, coeffx, x, y, z);
    const amrex::Real vx = dphidx_eval(degree, coeffy, x, y, z);
    const amrex::Real wx = dphidx_eval(degree, coeffz, x, y, z);

    return (ux * ux + uy * vx + uz * wx) - (ux * ux + uy * uy + uz * uz) +
           (3.0_rt * (ux * ux + vx * vx + wx * wx)) +
           (ux * ux + vx * uy + wx * uz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real nonlinearsum_12(
    int degree,
    const amrex::Real* coeffx,
    const amrex::Real* coeffy,
    const amrex::Real* coeffz,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    const amrex::Real ux = dphidx_eval(degree, coeffx, x, y, z);
    const amrex::Real uy = dphidy_eval(degree, coeffx, x, y, z);
    const amrex::Real uz = dphidz_eval(degree, coeffx, x, y, z);
    const amrex::Real vx = dphidx_eval(degree, coeffy, x, y, z);
    const amrex::Real vy = dphidy_eval(degree, coeffy, x, y, z);
    const amrex::Real vz = dphidz_eval(degree, coeffy, x, y, z);
    const amrex::Real wx = dphidx_eval(degree, coeffz, x, y, z);
    const amrex::Real wy = dphidy_eval(degree, coeffz, x, y, z);

    return (ux * uy + uy * vy + uz * wy) - (ux * vx + uy * vy + uz * vz) +
           (3.0_rt * (ux * uy + vx * vy + wx * wy)) +
           (ux * vx + vx * vy + wx * vz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real nonlinearsum_13(
    int degree,
    const amrex::Real* coeffx,
    const amrex::Real* coeffy,
    const amrex::Real* coeffz,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    const amrex::Real ux = dphidx_eval(degree, coeffx, x, y, z);
    const amrex::Real uy = dphidy_eval(degree, coeffx, x, y, z);
    const amrex::Real uz = dphidz_eval(degree, coeffx, x, y, z);
    const amrex::Real vx = dphidx_eval(degree, coeffy, x, y, z);
    const amrex::Real vz = dphidz_eval(degree, coeffy, x, y, z);
    const amrex::Real wx = dphidx_eval(degree, coeffz, x, y, z);
    const amrex::Real wy = dphidy_eval(degree, coeffz, x, y, z);
    const amrex::Real wz = dphidz_eval(degree, coeffz, x, y, z);

    return (ux * uz + uy * vz + uz * wz) - (ux * wx + uy * wy + uz * wz) +
           (3.0_rt * (ux * uz + vx * vz + wx * wz)) +
           (ux * wx + vx * wy + wx * wz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real nonlinearsum_21(
    int degree,
    const amrex::Real* coeffx,
    const amrex::Real* coeffy,
    const amrex::Real* coeffz,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    const amrex::Real ux = dphidx_eval(degree, coeffx, x, y, z);
    const amrex::Real uy = dphidy_eval(degree, coeffx, x, y, z);
    const amrex::Real uz = dphidz_eval(degree, coeffx, x, y, z);
    const amrex::Real vx = dphidx_eval(degree, coeffy, x, y, z);
    const amrex::Real vy = dphidy_eval(degree, coeffy, x, y, z);
    const amrex::Real vz = dphidz_eval(degree, coeffy, x, y, z);
    const amrex::Real wx = dphidx_eval(degree, coeffz, x, y, z);
    const amrex::Real wy = dphidy_eval(degree, coeffz, x, y, z);

    return (vx * ux + vy * vx + vz * wx) - (vx * ux + vy * uy + vz * uz) +
           (3.0_rt * (uy * ux + vy * vx + wy * wx)) +
           (uy * ux + vy * uy + wy * uz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real nonlinearsum_22(
    int degree,
    const amrex::Real* coeffx,
    const amrex::Real* coeffy,
    const amrex::Real* coeffz,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    const amrex::Real uy = dphidy_eval(degree, coeffx, x, y, z);
    const amrex::Real vx = dphidx_eval(degree, coeffy, x, y, z);
    const amrex::Real vy = dphidy_eval(degree, coeffy, x, y, z);
    const amrex::Real vz = dphidz_eval(degree, coeffy, x, y, z);
    const amrex::Real wy = dphidy_eval(degree, coeffz, x, y, z);

    return (vx * uy + vy * vy + vz * wy) - (vx * vx + vy * vy + vz * vz) +
           (3.0_rt * (uy * vx + vy * vy + wy * wy)) +
           (uy * vx + vy * vy + wy * vz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real nonlinearsum_23(
    int degree,
    const amrex::Real* coeffx,
    const amrex::Real* coeffy,
    const amrex::Real* coeffz,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    const amrex::Real uy = dphidy_eval(degree, coeffx, x, y, z);
    const amrex::Real uz = dphidz_eval(degree, coeffx, x, y, z);
    const amrex::Real vx = dphidx_eval(degree, coeffy, x, y, z);
    const amrex::Real vy = dphidy_eval(degree, coeffy, x, y, z);
    const amrex::Real vz = dphidz_eval(degree, coeffy, x, y, z);
    const amrex::Real wx = dphidx_eval(degree, coeffz, x, y, z);
    const amrex::Real wy = dphidy_eval(degree, coeffz, x, y, z);
    const amrex::Real wz = dphidz_eval(degree, coeffz, x, y, z);

    return (vx * uz + vy * vz + vz * wz) - (vx * wx + vy * wy + vz * wz) +
           (3.0_rt * (uy * wx + vy * vz + wy * wz)) +
           (uy * wx + vy * wy + wy * wz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real nonlinearsum_31(
    int degree,
    const amrex::Real* coeffx,
    const amrex::Real* coeffy,
    const amrex::Real* coeffz,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    const amrex::Real ux = dphidx_eval(degree, coeffx, x, y, z);
    const amrex::Real uy = dphidy_eval(degree, coeffx, x, y, z);
    const amrex::Real uz = dphidz_eval(degree, coeffx, x, y, z);
    const amrex::Real vx = dphidx_eval(degree, coeffy, x, y, z);
    const amrex::Real vz = dphidz_eval(degree, coeffy, x, y, z);
    const amrex::Real wx = dphidx_eval(degree, coeffz, x, y, z);
    const amrex::Real wy = dphidy_eval(degree, coeffz, x, y, z);
    const amrex::Real wz = dphidz_eval(degree, coeffz, x, y, z);

    return (wx * ux + wy * vx + wz * wx) - (wx * ux + wy * uy + wz * uz) +
           (3.0_rt * (uz * ux + vz * vx + wz * wx)) +
           (uz * ux + vz * uy + wz * uz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real nonlinearsum_32(
    int degree,
    const amrex::Real* coeffx,
    const amrex::Real* coeffy,
    const amrex::Real* coeffz,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    const amrex::Real uy = dphidy_eval(degree, coeffx, x, y, z);
    const amrex::Real uz = dphidz_eval(degree, coeffx, x, y, z);
    const amrex::Real vx = dphidx_eval(degree, coeffy, x, y, z);
    const amrex::Real vy = dphidy_eval(degree, coeffy, x, y, z);
    const amrex::Real vz = dphidz_eval(degree, coeffy, x, y, z);
    const amrex::Real wx = dphidx_eval(degree, coeffz, x, y, z);
    const amrex::Real wy = dphidy_eval(degree, coeffz, x, y, z);
    const amrex::Real wz = dphidz_eval(degree, coeffz, x, y, z);

    return (wx * uy + wy * vy + wz * wy) - (wx * vx + wy * vy + wz * vz) +
           (3.0_rt * (uz * uy + vz * vy + wz * wy)) +
           (uz * vx + vz * vy + wz * vz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real nonlinearsum_33(
    int degree,
    const amrex::Real* coeffx,
    const amrex::Real* coeffy,
    const amrex::Real* coeffz,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    const amrex::Real uz = dphidz_eval(degree, coeffx, x, y, z);
    const amrex::Real vz = dphidz_eval(degree, coeffy, x, y, z);
    const amrex::Real wx = dphidx_eval(degree, coeffz, x, y, z);
    const amrex::Real wy = dphidy_eval(degree, coeffz, x, y, z);
    const amrex::Real wz = dphidz_eval(degree, coeffz, x, y, z);

    return (wx * uz + wy * vz + wz * wz) - (wx * wx + wy * wy + wz * wz) +
           (3.0_rt * (uz * uz + vz * vz + wz * wz)) +
           (uz * wx + vz * wy + wz * wz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real vorticity_x(
    int degree,
    const amrex::Real* /* coeffx */,
    const amrex::Real* coeffy,
    const amrex::Real* coeffz,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    const amrex::Real vz = dphidz_eval(degree, coeffy, x, y, z);
    const amrex::Real wy = dphidy_eval(degree, coeffz, x, y, z);

    return wy - vz;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real vorticity_y(
    int degree,
    const amrex::Real* coeffx,
    const amrex::Real* /* coeffy */,
    const amrex::Real* coeffz,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    const amrex::Real uz = dphidz_eval(degree, coeffx, x, y, z);
    const amrex::Real wx = dphidx_eval(degree, coeffz, x, y, z);

    return uz - wx;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real vorticity_z(
    int degree,
    const amrex::Real* coeffx,
    const amrex::Real* coeffy,
    const amrex::Real* /* coeffz */,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    const amrex::Real uy = dphidy_eval(degree, coeffx, x, y, z);
    const amrex::Real vx = dphidx_eval(degree, coeffy, x, y, z);

    return vx - uy;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real vorticity_mag(
    int degree,
    const amrex::Real* coeffx,
    const amrex::Real* coeffy,
    const amrex::Real* coeffz,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    const amrex::Real uy = dphidy_eval(degree, coeffx, x, y, z);
    const amrex::Real uz = dphidz_eval(degree, coeffx, x, y, z);

    const amrex::Real vx = dphidx_eval(degree, coeffy, x, y, z);
    const amrex::Real vz = dphidz_eval(degree, coeffy, x, y, z);

    const amrex::Real wx = dphidx_eval(degree, coeffz, x, y, z);
    const amrex::Real wy = dphidy_eval(degree, coeffz, x, y, z);

    return std::sqrt(
        ((uy - vx) * (uy - vx)) + ((vz - wy) * (vz - wy)) +
        ((wx - uz) * (wx - uz)));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real curvature(
    int degree,
    const amrex::Real* coeff,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    const amrex::Real phix = dphidx_eval(degree, coeff, x, y, z);
    const amrex::Real phiy = dphidy_eval(degree, coeff, x, y, z);
    const amrex::Real phiz = dphidz_eval(degree, coeff, x, y, z);

    amrex::Real phixx = 0.0_rt;
    amrex::Real phixy = 0.0_rt;
    amrex::Real phiyy = 0.0_rt;
    amrex::Real phixz = 0.0_rt;
    amrex::Real phiyz = 0.0_rt;
    amrex::Real phizz = 0.0_rt;

    int ind = 0;
    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j) {
            for (int k = 0; k <= degree; ++k) {

                phixx += i * amrex::max(i - 1, 0) * coeff[ind] *
                         static_cast<amrex::Real>(
                             std::pow(x, amrex::max(i - 2, 0))) *
                         static_cast<amrex::Real>(std::pow(y, j)) *
                         static_cast<amrex::Real>(std::pow(z, k));

                phixy += i * j * coeff[ind] *
                         static_cast<amrex::Real>(
                             std::pow(x, amrex::max(i - 1, 0))) *
                         static_cast<amrex::Real>(
                             std::pow(y, amrex::max(j - 1, 0))) *
                         static_cast<amrex::Real>(std::pow(z, k));

                phiyy += j * amrex::max(j - 1, 0) * coeff[ind] *
                         static_cast<amrex::Real>(std::pow(x, i)) *
                         static_cast<amrex::Real>(
                             std::pow(y, amrex::max(j - 2, 0))) *
                         static_cast<amrex::Real>(std::pow(z, k));

                phixz +=
                    i * k * coeff[ind] *
                    static_cast<amrex::Real>(
                        std::pow(x, amrex::max(i - 1, 0))) *
                    static_cast<amrex::Real>(std::pow(y, j)) *
                    static_cast<amrex::Real>(std::pow(z, amrex::max(k - 1, 0)));

                phiyz +=
                    j * k * coeff[ind] *
                    static_cast<amrex::Real>(std::pow(x, i)) *
                    static_cast<amrex::Real>(
                        std::pow(y, amrex::max(j - 1, 0))) *
                    static_cast<amrex::Real>(std::pow(z, amrex::max(k - 1, 0)));

                phizz +=
                    k * amrex::max(k - 1, 0) * coeff[ind] *
                    static_cast<amrex::Real>(std::pow(x, i)) *
                    static_cast<amrex::Real>(std::pow(y, j)) *
                    static_cast<amrex::Real>(std::pow(z, amrex::max(k - 2, 0)));

                ++ind;
            }
        }
    }

    return -((phix * phix * phiyy) - (2.0_rt * phix * phiy * phixy) +
             (phiy * phiy * phixx) + (phix * phix * phizz) -
             (2.0_rt * phix * phiz * phixz) + (phiz * phiz * phixx) +
             (phiy * phiy * phizz) - (2.0_rt * phiy * phiz * phiyz) +
             (phiz * phiz * phiyy)) /
           std::pow((phix * phix) + (phiy * phiy) + (phiz * phiz), 1.5_rt);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real q_criterion(
    int degree,
    const amrex::Real* coeffx,
    const amrex::Real* coeffy,
    const amrex::Real* coeffz,
    amrex::Real x,
    amrex::Real y,
    amrex::Real z)
{
    const amrex::Real ux = dphidx_eval(degree, coeffx, x, y, z);
    const amrex::Real uy = dphidy_eval(degree, coeffx, x, y, z);
    const amrex::Real uz = dphidz_eval(degree, coeffx, x, y, z);

    const amrex::Real vx = dphidx_eval(degree, coeffy, x, y, z);
    const amrex::Real vy = dphidy_eval(degree, coeffy, x, y, z);
    const amrex::Real vz = dphidz_eval(degree, coeffy, x, y, z);

    const amrex::Real wx = dphidx_eval(degree, coeffz, x, y, z);
    const amrex::Real wy = dphidy_eval(degree, coeffz, x, y, z);
    const amrex::Real wz = dphidz_eval(degree, coeffz, x, y, z);

    const amrex::Real S2 = (2.0_rt * ux * ux) + (2.0_rt * vy * vy) +
                           (2.0_rt * wz * wz) + ((uy + vx) * (uy + vx)) +
                           ((vz + wy) * (vz + wy)) + ((wx + uz) * (wx + uz));
    const amrex::Real W2 = ((uy - vx) * (uy - vx)) + ((vz - wy) * (vz - wy)) +
                           ((wx - uz) * (wx - uz));

    return 0.25_rt * (W2 - S2);
}

} // namespace amr_wind_tests::analytical_function
#endif
