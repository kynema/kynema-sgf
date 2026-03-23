#ifndef SLICE_H
#define SLICE_H

#include <iterator>
#include "AMReX_Vector.H"
#include "AMReX_Gpu.H"

namespace amr_wind::utils {

template <typename T>
struct Slice
{
    using value_type = T;
    using size_t = std::size_t;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reference = T&;
    using pointer = T*;
    using iterator = T*;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reference = const T&;
    using const_pointer = const T*;
    using const_iterator = const T*;
    using const_reverse_iterator = const std::reverse_iterator<iterator>;

    pointer ptr_begin;
    pointer ptr_end;

    explicit Slice() : ptr_begin(nullptr), ptr_end(nullptr) {}

    explicit Slice(T* pbegin, const size_t n)
        : ptr_begin(pbegin), ptr_end(pbegin + n)
    {}

    explicit Slice(T* pbegin, T* pend) : ptr_begin(pbegin), ptr_end(pend) {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE reference
    operator[](const size_t idx)
    {
        return *(ptr_begin + idx);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE const_reference
    operator[](const size_type idx) const
    {
        return *(ptr_begin + idx);
    }

    // NOLINTNEXTLINE(modernize-use-nodiscard)
    AMREX_GPU_HOST_DEVICE size_type size() const
    {
        return (ptr_end - ptr_begin);
    }

    AMREX_GPU_HOST_DEVICE pointer data() { return ptr_begin; }
    AMREX_GPU_HOST_DEVICE iterator begin() { return ptr_begin; }
    AMREX_GPU_HOST_DEVICE iterator end() { return ptr_end; }
    // NOLINTNEXTLINE(modernize-use-nodiscard)
    AMREX_GPU_HOST_DEVICE const_iterator data() const { return ptr_begin; }
    // NOLINTNEXTLINE(modernize-use-nodiscard)
    AMREX_GPU_HOST_DEVICE const_iterator begin() const { return ptr_begin; }
    // NOLINTNEXTLINE(modernize-use-nodiscard)
    AMREX_GPU_HOST_DEVICE const_iterator end() const { return ptr_end; }
    // NOLINTNEXTLINE(modernize-use-nodiscard)
    AMREX_GPU_HOST_DEVICE const_iterator cbegin() const { return ptr_begin; }
    // NOLINTNEXTLINE(modernize-use-nodiscard)
    AMREX_GPU_HOST_DEVICE const_iterator cend() const { return ptr_end; }
};

template <typename T>
Slice<T> slice(std::vector<T>& vec, const size_t start, const size_t count)
{
    AMREX_ASSERT((start + count) <= vec.size());
    return Slice<T>{&vec[start], count};
}

template <typename T>
Slice<T> slice(std::vector<T>& vec, const size_t start)
{
    AMREX_ASSERT(start < vec.size());
    return slice(vec, start, vec.size() - start);
}

template <typename T>
Slice<T> slice(amrex::Vector<T>& vec, const size_t start, const size_t count)
{
    AMREX_ASSERT((start + count) <= static_cast<size_t>(vec.size()));
    return Slice<T>{&vec[start], count};
}

template <typename T>
Slice<T> slice(amrex::Vector<T>& vec, const size_t start)
{
    AMREX_ASSERT(start < static_cast<size_t>(vec.size()));
    return slice(vec, start, vec.size() - start);
}

template <typename T>
Slice<T>
slice(amrex::Gpu::DeviceVector<T>& vec, const size_t start, const size_t count)
{
    AMREX_ASSERT((start + count) <= static_cast<size_t>(vec.size()));
    return Slice<T>{(vec.data() + start), count};
}

template <typename T>
Slice<T> slice(amrex::Gpu::DeviceVector<T>& vec, const size_t start)
{
    AMREX_ASSERT(start < static_cast<size_t>(vec.size()));
    return slice(vec, start, vec.size() - start);
}

} // namespace amr_wind::utils

#endif /* SLICE_H */
