#ifndef VS_VECTORI_H
#define VS_VECTORI_H

#include <ostream>
#include <cmath>
#include "amr-wind/core/vs/vector.H"

namespace amr_wind::vs {

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T>
VectorT<T>::operator-() const
{
    return VectorT<T>{-vv[0], -vv[1], -vv[2]};
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T>
VectorT<T>::operator*=(const T fac)
{
    vv[0] *= fac;
    vv[1] *= fac;
    vv[2] *= fac;
    return *this;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T>
VectorT<T>::operator*=(const VectorT<T>& vin)
{
    vv[0] *= vin.x();
    vv[1] *= vin.y();
    vv[2] *= vin.z();
    return *this;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T>
VectorT<T>::operator/=(const T fac)
{
    vv[0] /= fac;
    vv[1] /= fac;
    vv[2] /= fac;
    return *this;
}

template <typename T, typename OStream>
OStream& operator<<(OStream& out, const VectorT<T>& vec)
{
    out << "(" << vec.x() << " " << vec.y() << " " << vec.z() << ")";
    return out;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T>
operator+(const VectorT<T>& v1, const VectorT<T>& v2)
{
    return VectorT<T>{v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z()};
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T>
operator-(const VectorT<T>& v1, const VectorT<T>& v2)
{
    return VectorT<T>{v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z()};
}

template <typename T1, typename T2>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T1>
operator*(const VectorT<T1>& inp, const T2 fac)
{
    return VectorT<T1>{inp.x() * fac, inp.y() * fac, inp.z() * fac};
}

template <typename T1, typename T2>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T1>
operator*(const T2 fac, const VectorT<T1>& inp)
{
    return VectorT<T1>{inp.x() * fac, inp.y() * fac, inp.z() * fac};
}

template <typename T1, typename T2>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T1>
operator/(const VectorT<T1>& inp, const T2 fac)
{
    return VectorT<T1>{inp.x() / fac, inp.y() / fac, inp.z() / fac};
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T
operator&(const VectorT<T>& v1, const VectorT<T>& v2)
{
    return (v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z());
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T>
operator^(const VectorT<T>& v1, const VectorT<T>& v2)
{
    return VectorT<T>{
        (v1.y() * v2.z() - v1.z() * v2.y()),
        (v1.z() * v2.x() - v1.x() * v2.z()),
        (v1.x() * v2.y() - v1.y() * v2.x())};
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T>
operator*(const VectorT<T>& v1, const VectorT<T>& v2)
{
    return VectorT<T>{v1.x() * v2.x(), v1.y() * v2.y(), v1.z() * v2.z()};
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T mag_sqr(const VectorT<T>& v)
{
    return (v.x() * v.x() + v.y() * v.y() + v.z() * v.z());
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T mag(const VectorT<T>& v)
{
    return std::sqrt(mag_sqr(v));
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T>& VectorT<T>::normalize()
{
    T vmag = mag(*this);

    if (vmag < Traits::eps()) {
        *this = VectorT<T>::zero();
    } else {
        *this /= vmag;
    }
    return *this;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> sqrt(const VectorT<T>& v)
{
    return VectorT<T>{std::sqrt(v.x()), std::sqrt(v.y()), std::sqrt(v.z())};
}

} // namespace amr_wind::vs

#endif /* VS_VECTORI_H */
