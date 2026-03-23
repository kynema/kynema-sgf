#ifndef VS_TENSORI_H
#define VS_TENSORI_H

#include <cmath>
#include "amr-wind/core/vs/tensor.H"
#include "amr-wind/utilities/trig_ops.H"

namespace amr_wind::vs {

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE TensorT<T>::TensorT(
    const VectorT<T>& x,
    const VectorT<T>& y,
    const VectorT<T>& z,
    const bool transpose)
{
    if (transpose) {
        this->cols(x, y, z);
    } else {
        this->rows(x, y, z);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
TensorT<T>::rows(const VectorT<T>& x, const VectorT<T>& y, const VectorT<T>& z)
{
    // clang-format off
    vv[0] = x.x(); vv[1] = x.y(); vv[2] = x.z();
    vv[3] = y.x(); vv[4] = y.y(); vv[5] = y.z();
    vv[6] = z.x(); vv[7] = z.y(); vv[8] = z.z();
    // clang-format on
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
TensorT<T>::cols(const VectorT<T>& x, const VectorT<T>& y, const VectorT<T>& z)
{
    // clang-format off
    vv[0] = x.x(); vv[1] = y.x(); vv[2] = z.x();
    vv[3] = x.y(); vv[4] = y.y(); vv[5] = z.y();
    vv[6] = x.z(); vv[7] = y.z(); vv[8] = z.z();
    // clang-format on
}

// clang-format off
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
VectorT<T> TensorT<T>::x() const 
{
    return VectorT<T>{vv[0], vv[1], vv[2]};
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
VectorT<T> TensorT<T>::y() const 
{
    return VectorT<T>{vv[3], vv[4], vv[5]};
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
VectorT<T> TensorT<T>::z() const 
{
    return VectorT<T>{vv[6], vv[7], vv[8]};
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
VectorT<T> TensorT<T>::cx() const 
{
    return VectorT<T>{vv[0], vv[3], vv[6]};
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
VectorT<T> TensorT<T>::cy() const 
{
    return VectorT<T>{vv[1], vv[4], vv[7]};
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
VectorT<T> TensorT<T>::cz() const 
{
    return VectorT<T>{vv[2], vv[5], vv[8]};
}
// clang-format on

template <typename T, typename OStream>
OStream& operator<<(OStream& out, const TensorT<T>& t)
{
    out << "(";
    for (int i = 0; i < t.ncomp; ++i) {
        out << " " << t[i];
    }
    out << " )";
    return out;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T>
operator&(const TensorT<T>& t, const VectorT<T>& v)
{
    VectorT<T> out;
    out.x() = t.xx() * v.x() + t.xy() * v.y() + t.xz() * v.z();
    out.y() = t.yx() * v.x() + t.yy() * v.y() + t.yz() * v.z();
    out.z() = t.zx() * v.x() + t.zy() * v.y() + t.zz() * v.z();
    return out;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T>
operator&(const VectorT<T>& v, const TensorT<T>& t)
{
    VectorT<T> out;
    out.x() = t.xx() * v.x() + t.yx() * v.y() + t.zx() * v.z();
    out.y() = t.xy() * v.x() + t.yy() * v.y() + t.zy() * v.z();
    out.z() = t.xz() * v.x() + t.yz() * v.y() + t.zz() * v.z();
    return out;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE TensorT<T>
operator&(const TensorT<T>& t1, const TensorT<T>& t2)
{
    TensorT<T> t3;

    t3.xx() = t1.xx() * t2.xx() + t1.xy() * t2.yx() + t1.xz() * t2.zx();
    t3.xy() = t1.xx() * t2.xy() + t1.xy() * t2.yy() + t1.xz() * t2.zy();
    t3.xz() = t1.xx() * t2.xz() + t1.xy() * t2.yz() + t1.xz() * t2.zz();

    t3.yx() = t1.yx() * t2.xx() + t1.yy() * t2.yx() + t1.yz() * t2.zx();
    t3.yy() = t1.yx() * t2.xy() + t1.yy() * t2.yy() + t1.yz() * t2.zy();
    t3.yz() = t1.yx() * t2.xz() + t1.yy() * t2.yz() + t1.yz() * t2.zz();

    t3.zx() = t1.zx() * t2.xx() + t1.zy() * t2.yx() + t1.zz() * t2.zx();
    t3.zy() = t1.zx() * t2.xy() + t1.zy() * t2.yy() + t1.zz() * t2.zy();
    t3.zz() = t1.zx() * t2.xz() + t1.zy() * t2.yz() + t1.zz() * t2.zz();

    return t3;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE TensorT<T>
operator+(const TensorT<T>& t1, const TensorT<T>& t2)
{
    return TensorT<T>{
        t1.vv[0] + t2.vv[0], t1.vv[1] + t2.vv[1], t1.vv[2] + t2.vv[2],
        t1.vv[3] + t2.vv[3], t1.vv[4] + t2.vv[4], t1.vv[5] + t2.vv[5],
        t1.vv[6] + t2.vv[6], t1.vv[7] + t2.vv[7], t1.vv[8] + t2.vv[8]};
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE TensorT<T>
operator-(const TensorT<T>& t1, const TensorT<T>& t2)
{
    return TensorT<T>{
        t1.vv[0] - t2.vv[0], t1.vv[1] - t2.vv[1], t1.vv[2] - t2.vv[2],
        t1.vv[3] - t2.vv[3], t1.vv[4] - t2.vv[4], t1.vv[5] - t2.vv[5],
        t1.vv[6] - t2.vv[6], t1.vv[7] - t2.vv[7], t1.vv[8] - t2.vv[8]};
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T
operator&&(const TensorT<T>& t1, const TensorT<T>& t2)
{
    return (
        t1.vv[0] * t2.vv[0] + t1.vv[1] * t2.vv[1] + t1.vv[2] * t2.vv[2] +
        t1.vv[3] * t2.vv[3] + t1.vv[4] * t2.vv[4] + t1.vv[5] * t2.vv[5] +
        t1.vv[6] * t2.vv[6] + t1.vv[7] * t2.vv[7] + t1.vv[8] * t2.vv[8]);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T mag_sqr(const TensorT<T>& t)
{
    // cppcheck-suppress duplicateExpression
    return (t && t);
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T mag(const TensorT<T>& t)
{
    return std::sqrt(mag_sqr(t));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor xrot(const amrex::Real angle)
{
    const amrex::Real ang = utils::radians(angle);
    const amrex::Real cval = std::cos(ang);
    const amrex::Real sval = std::sin(ang);

    return Tensor{1.0_rt, 0.0_rt, 0.0_rt, 0.0_rt, cval,
                  sval,   0.0_rt, -sval,  cval};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor yrot(const amrex::Real angle)
{
    const amrex::Real ang = utils::radians(angle);
    const amrex::Real cval = std::cos(ang);
    const amrex::Real sval = std::sin(ang);

    return Tensor{cval,   0.0_rt, -sval,  0.0_rt, 1.0_rt,
                  0.0_rt, sval,   0.0_rt, cval};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor zrot(const amrex::Real angle)
{
    const amrex::Real ang = utils::radians(angle);
    const amrex::Real cval = std::cos(ang);
    const amrex::Real sval = std::sin(ang);

    return Tensor{cval,   sval,   0.0_rt, -sval, cval,
                  0.0_rt, 0.0_rt, 0.0_rt, 1.0_rt};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor
quaternion(const Vector& axis, const amrex::Real angle)
{
    const amrex::Real ang = -1.0_rt * utils::radians(angle);
    const amrex::Real cval = std::cos(0.5_rt * ang);
    const amrex::Real sval = std::sin(0.5_rt * ang);
    const auto vmag = mag(axis);
    const amrex::Real q0 = cval;
    const amrex::Real q1 = sval * axis.x() / vmag;
    const amrex::Real q2 = sval * axis.y() / vmag;
    const amrex::Real q3 = sval * axis.z() / vmag;

    Tensor t;
    t.xx() = (q0 * q0) + (q1 * q1) - (q2 * q2) - (q3 * q3);
    t.xy() = 2.0_rt * (q1 * q2 - q0 * q3);
    t.xz() = 2.0_rt * (q0 * q2 + q1 * q3);

    t.yx() = 2.0_rt * (q1 * q2 + q0 * q3);
    t.yy() = (q0 * q0) - (q1 * q1) + (q2 * q2) - (q3 * q3);
    t.yz() = 2.0_rt * (q2 * q3 - q0 * q1);

    t.zx() = 2.0_rt * (q1 * q3 - q0 * q2);
    t.zy() = 2.0_rt * (q0 * q1 + q2 * q3);
    t.zz() = (q0 * q0) - (q1 * q1) - (q2 * q2) + (q3 * q3);

    return t;
}

} // namespace amr_wind::vs

#endif /* VS_TENSORI_H */
