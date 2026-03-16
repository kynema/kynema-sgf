#ifndef TENSOR_OPS_H
#define TENSOR_OPS_H

#include "AMReX_Gpu.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::utils {

/** Compute the dot product of two vectors
 */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T dot_prod(const T* x, const T* y)
{
    static_assert(
        AMREX_SPACEDIM == 3, "Vector ops only defined for AMREX_SPACEDIM == 3");
    return ((x[0] * y[0]) + (x[1] * y[1]) + (x[2] * y[2]));
}

/** Return the magnitude of a 3-D vector
 */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T vec_mag(const T* x)
{
    return std::sqrt(dot_prod(x, x));
}

/** Convert a given vector to a unit vector
 */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void vec_normalize(T* x)
{
    static_assert(
        AMREX_SPACEDIM == 3, "Vector ops only defined for AMREX_SPACEDIM == 3");
    const T inv_mag = 1.0_rt / vec_mag(x);
    x[0] *= inv_mag;
    x[1] *= inv_mag;
    x[2] *= inv_mag;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T
dot_prod(const T* x, const T* y, const int ndim)
{
    T sum = 0.0_rt;
    for (int i = 0; i < ndim; ++i) {
        sum += x[i] * y[i];
    }
    return sum;
}

/** Compute the cross product of two vectors `a`, `b` and update `c`
 */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cross_prod(const T* a, const T* b, T* c)
{
    static_assert(
        AMREX_SPACEDIM == 3, "Vector ops only defined for AMREX_SPACEDIM == 3");
    c[0] = (a[1] * b[2]) - (a[2] * b[1]);
    c[1] = (a[2] * b[0]) - (a[0] * b[2]);
    c[2] = (a[0] * b[1]) - (a[1] * b[0]);
}

/** Compute the transformation `{y} = [A] {x}`
 */
template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void transform_vec(
    const amrex::Array2D<T, 0, AMREX_SPACEDIM, 0, AMREX_SPACEDIM>& tmat,
    const T* x,
    T* y)
{
    static_assert(
        AMREX_SPACEDIM == 3, "Vector ops only defined for AMREX_SPACEDIM == 3");
    y[0] = (tmat(0, 0) * x[0]) + (tmat(0, 1) * x[1]) + (tmat(0, 2) * x[2]);
    y[1] = (tmat(1, 0) * x[0]) + (tmat(1, 1) * x[1]) + (tmat(1, 2) * x[2]);
    y[2] = (tmat(2, 0) * x[0]) + (tmat(2, 1) * x[1]) + (tmat(2, 2) * x[2]);
}

/** Compute the transformation `{y} = {x}^T [A]`
 */
template <typename T>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void inv_transform_vec(
    const amrex::Array2D<T, 0, AMREX_SPACEDIM, 0, AMREX_SPACEDIM>& tmat,
    const T* x,
    T* y)
{
    static_assert(
        AMREX_SPACEDIM == 3, "Vector ops only defined for AMREX_SPACEDIM == 3");
    y[0] = (tmat(0, 0) * x[0]) + (tmat(1, 0) * x[1]) + (tmat(2, 0) * x[2]);
    y[1] = (tmat(0, 1) * x[0]) + (tmat(1, 1) * x[1]) + (tmat(2, 1) * x[2]);
    y[2] = (tmat(0, 2) * x[0]) + (tmat(1, 2) * x[1]) + (tmat(2, 2) * x[2]);
}

} // namespace amr_wind::utils

#endif /* TENSOR_OPS_H */
