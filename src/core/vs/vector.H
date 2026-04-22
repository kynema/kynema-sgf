#ifndef VS_VECTOR_H
#define VS_VECTOR_H

#include "AMReX_Gpu.H"
#include "amr-wind/core/vs/vstraits.H"

namespace amr_wind::vs {

/** Vector in 3D space
 */
template <typename T>
struct VectorT
{
    amrex::GpuArray<T, 3> vv = {Traits::zero(), Traits::zero(), Traits::zero()};

    //! Number of components
    static constexpr int ncomp = 3;
    using size_type = int;
    using value_type = T;
    using reference = T&;
    using iterator = T*;
    using const_iterator = const T*;
    using Traits = DTraits<T>;
    using VType = VectorT<T>;

    //! Construct a default vector, all components set to zero
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT() = default;

    /** New vector given the three components
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    VectorT(const T& x, const T& y, const T& z)
        : vv{x, y, z}
    {}

    ~VectorT() = default;
    VectorT(const VectorT&) = default;
    VectorT(VectorT&&) = default;
    VectorT& operator=(const VectorT&) & = default;
    VectorT& operator=(const VectorT&) && = delete;
    VectorT& operator=(VectorT&&) & = default;
    VectorT& operator=(VectorT&&) && = delete;

    //! Zero vector
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static constexpr VectorT<T> zero()
    {
        return VectorT<T>{Traits::zero(), Traits::zero(), Traits::zero()};
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static constexpr VectorT<T> one()
    {
        return VectorT<T>{Traits::one(), Traits::one(), Traits::one()};
    }

    /** Vector along x-axis
     *
     *  \param x Magnitude of vector
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static constexpr VectorT<T>
    ihat(const T& x = Traits::one())
    {
        return VectorT<T>{x, Traits::zero(), Traits::zero()};
    }

    /** Vector along y-axis
     *
     *  \param y Magnitude of vector
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static constexpr VectorT<T>
    jhat(const T& y = Traits::one())
    {
        return VectorT<T>{Traits::zero(), y, Traits::zero()};
    }

    /** Vector along z-axis
     *
     *  \param z Magnitude of vector
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static constexpr VectorT<T>
    khat(const T& z = Traits::one())
    {
        return VectorT<T>{Traits::zero(), Traits::zero(), z};
    }

    //! Normalize the vector to unit vector
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T>& normalize();

    //! Sqrt of each component
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> sqrt(const VectorT&);

    //! Return the unit vector parallel to this vector
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> unit() const
    {
        return VectorT<T>(*this).normalize();
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& x() & { return vv[0]; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& y() & { return vv[1]; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& z() & { return vv[2]; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& x() const&
    {
        return vv[0];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& y() const&
    {
        return vv[1];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& z() const&
    {
        return vv[2];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> operator-() const;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> operator*=(T val);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T>
    operator*=(const VectorT<T>& valvec);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> operator/=(T val);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& operator[](size_type pos) &
    {
        return vv[pos];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T&
    operator[](size_type pos) const&
    {
        return vv[pos];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T* data() { return vv.data(); }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T* data() const
    {
        return vv.data();
    }

    iterator begin() { return vv.begin(); }
    iterator end() { return vv.end(); }
    [[nodiscard]] const_iterator cbegin() const { return vv.cbegin(); }
    [[nodiscard]] const_iterator cend() const { return vv.cend(); }
    [[nodiscard]] size_type size() const { return ncomp; }
};

using Vector = VectorT<amrex::Real>;

} // namespace amr_wind::vs

#include "amr-wind/core/vs/vectorI.H"

#endif /* VS_VECTOR_H */
