#ifndef VS_TENSOR_H
#define VS_TENSOR_H

#include "AMReX_Gpu.H"
#include "amr-wind/core/vs/vstraits.H"
#include "amr-wind/core/vs/vector.H"

namespace amr_wind::vs {

/** Tensor in 3D vector space
 */
template <typename T>
struct TensorT
{
    amrex::GpuArray<T, 9> vv = {Traits::zero(), Traits::zero(), Traits::zero(),
                                Traits::zero(), Traits::zero(), Traits::zero(),
                                Traits::zero(), Traits::zero(), Traits::zero()};

    static constexpr int ncomp = 9;
    using size_type = int;
    using value_type = T;
    using reference = T&;
    using iterator = T*;
    using const_iterator = const T*;
    using Traits = DTraits<T>;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr TensorT() = default;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr TensorT(
        const T& xx,
        const T& xy,
        const T& xz,
        const T& yx,
        const T& yy,
        const T& yz,
        const T& zx,
        const T& zy,
        const T& zz)
        : vv{xx, xy, xz, yx, yy, yz, zx, zy, zz}
    {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE TensorT(
        const VectorT<T>& x,
        const VectorT<T>& y,
        const VectorT<T>& z,
        bool transpose = false);

    ~TensorT() = default;
    TensorT(const TensorT&) = default;
    TensorT(TensorT&&) = default;
    TensorT& operator=(const TensorT&) & = default;
    TensorT& operator=(const TensorT&) && = delete;
    TensorT& operator=(TensorT&&) & = default;
    TensorT& operator=(TensorT&&) && = delete;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static constexpr TensorT<T> zero()
    {
        return TensorT<T>{Traits::zero(), Traits::zero(), Traits::zero(),
                          Traits::zero(), Traits::zero(), Traits::zero(),
                          Traits::zero(), Traits::zero(), Traits::zero()};
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static constexpr TensorT<T>
    identity()
    {
        return TensorT{Traits::one(),  Traits::zero(), Traits::zero(),
                       Traits::zero(), Traits::one(),  Traits::zero(),
                       Traits::zero(), Traits::zero(), Traits::one()};
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    rows(const VectorT<T>& x, const VectorT<T>& y, const VectorT<T>& z);
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
    cols(const VectorT<T>& x, const VectorT<T>& y, const VectorT<T>& z);

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> x() const;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> y() const;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> z() const;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> cx() const;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> cy() const;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE VectorT<T> cz() const;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& xx() & { return vv[0]; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& xy() & { return vv[1]; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& xz() & { return vv[2]; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& yx() & { return vv[3]; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& yy() & { return vv[4]; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& yz() & { return vv[5]; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& zx() & { return vv[6]; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& zy() & { return vv[7]; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE T& zz() & { return vv[8]; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& xx() const&
    {
        return vv[0];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& xy() const&
    {
        return vv[1];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& xz() const&
    {
        return vv[2];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& yx() const&
    {
        return vv[3];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& yy() const&
    {
        return vv[4];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& yz() const&
    {
        return vv[5];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& zx() const&
    {
        return vv[6];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& zy() const&
    {
        return vv[7];
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const T& zz() const&
    {
        return vv[8];
    }

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

using Tensor = TensorT<amrex::Real>;

} // namespace amr_wind::vs

#include "amr-wind/core/vs/tensorI.H"

#endif /* VS_TENSOR_H */
