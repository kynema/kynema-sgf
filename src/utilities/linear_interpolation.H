#ifndef LINEAR_INTERPOLATION_H
#define LINEAR_INTERPOLATION_H

#include "AMReX_Gpu.H"
#include "AMReX_REAL.H"
#include "AMReX_Extension.H"

using namespace amrex::literals;

namespace amr_wind::interp {

enum class Limits : int { // NOLINT(performance-enum-size)
    LOWLIM = -2,
    UPLIM = -1,
    VALID = 0,
};

struct Index
{
    int idx;
    Limits lim;
};

template <typename It, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Index
check_bounds(const It begin, const It end, const T& x)
{
    const int sz = static_cast<int>(end - begin);

    if ((sz < 2) || (x < *begin)) {
        return Index{.idx = 0, .lim = Limits::LOWLIM};
    }
    if (x > *(begin + (sz - 1))) {
        return Index{.idx = sz - 1, .lim = Limits::UPLIM};
    }

    return Index{.idx = 0, .lim = Limits::VALID};
}

template <typename It, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Index
bisection_search(const It begin, const It end, const T& x)
{
    auto idx = check_bounds(begin, end, x);
    if ((idx.lim == Limits::LOWLIM) || (idx.lim == Limits::UPLIM)) {
        return idx;
    }

    int il = 0;
    int ir = end - begin;
    const T xl = *begin;

    while ((ir - il) > 1) {
        int mid = (il + ir) >> 1;
        const T xmid = *(begin + mid);

        if ((x - xmid) * (x - xl) <= 0.0_rt) {
            ir = mid;
        } else {
            il = mid;
        }
    }
    idx.idx = il;
    return idx;
}

template <typename It, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Index
nearest_search(const It begin, const It end, const T& x)
{

    auto idx = bisection_search(begin, end, x);
    if ((idx.lim == Limits::LOWLIM) || (idx.lim == Limits::UPLIM)) {
        return idx;
    }

    if ((x - begin[idx.idx]) > (begin[idx.idx + 1] - x)) {
        idx.idx += 1;
    }

    return idx;
}

template <typename It, typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Index
find_index(const It begin, const It end, const T& x, const int hint = 1)
{
    auto idx = check_bounds(begin, end, x);
    if ((idx.lim == Limits::LOWLIM) || (idx.lim == Limits::UPLIM)) {
        return idx;
    }

    for (It it = (begin + hint); it < end; ++it) {
        if (x <= *it) {
            idx.idx = it - begin - 1;
            break;
        }
    }
    return idx;
}

template <typename C1, typename C2>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    typename std::iterator_traits<C2>::value_type
    linear_impl(
        const C1 xbegin,
        const C2 yinp,
        const typename std::iterator_traits<C1>::value_type& xout,
        const Index& idx,
        const int ncomp = 1,
        const int comp = 0)
{
    using DType1 = typename std::iterator_traits<C1>::value_type;

    if ((idx.lim == Limits::LOWLIM) || (idx.lim == Limits::UPLIM)) {
        return yinp[(ncomp * idx.idx) + comp];
    }
    static constexpr DType1 eps = std::numeric_limits<float>::epsilon();
    const int j = idx.idx;
    const auto denom = (xbegin[j + 1] - xbegin[j]);
    const auto facR = (denom > eps) ? ((xout - xbegin[j]) / denom) : 1.0_rt;
    const auto facL = static_cast<DType1>(1.0_rt) - facR;
    return (facL * yinp[(ncomp * j) + comp]) +
           (facR * yinp[(ncomp * (j + 1)) + comp]);
}

template <typename C1, typename C2>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    typename std::iterator_traits<C2>::value_type
    linear(
        const C1 xbegin,
        const C1 xend,
        const C2 yinp,
        const typename std::iterator_traits<C1>::value_type& xout,
        const int ncomp = 1,
        const int comp = 0)
{
    const auto idx = bisection_search(xbegin, xend, xout);
    return linear_impl(xbegin, yinp, xout, idx, ncomp, comp);
}

template <typename C1, typename C2>
typename C2::value_type linear(
    const C1& xinp,
    const C2& yinp,
    const typename C1::value_type& xout,
    const int ncomp = 1,
    const int comp = 0)
{
    return linear(
        xinp.data(), (xinp.data() + xinp.size()), yinp.data(), xout, ncomp,
        comp);
}

template <typename C1, typename C2>
void linear_monotonic(
    const C1& xinp,
    const C2& yinp,
    const C1& xout,
    C2& yout,
    const int ncomp = 1,
    const int comp = 0)
{
    static constexpr typename C1::value_type eps =
        std::numeric_limits<float>::epsilon();
    AMREX_ASSERT(xinp.size() == yinp.size());
    AMREX_ASSERT(xout.size() == yout.size());

    int hint = 1;
    int npts = xout.size();
    for (int i = 0; i < npts; ++i) {
        const auto& x = xout[i];
        const auto idx =
            find_index(xinp.data(), (xinp.data() + xinp.size()), x, hint);

        if ((idx.lim == Limits::LOWLIM) || (idx.lim == Limits::UPLIM)) {
            yout[i] = yinp[(ncomp * idx.idx) + comp];
        } else if (idx.lim == Limits::VALID) {
            int j = idx.idx;
            const auto denom = (xinp[j + 1] - xinp[j]);
            const auto facR = (denom > eps) ? ((x - xinp[j]) / denom) : 1.0_rt;
            const auto facL =
                static_cast<typename C1::value_type>(1.0_rt) - facR;
            yout[i] = (facL * yinp[(ncomp * j) + comp]) +
                      (facR * yinp[(ncomp * (j + 1)) + comp]);
            ;
        }
        hint = idx.idx + 1;
    }
}

template <typename C1, typename C2>
void linear(
    const C1& xinp,
    const C2& yinp,
    const C1& xout,
    C2& yout,
    const int ncomp = 1,
    const int comp = 0)
{
    AMREX_ASSERT(
        static_cast<amrex::Long>(xinp.size()) ==
        static_cast<amrex::Long>(yinp.size()));
    AMREX_ASSERT(
        static_cast<amrex::Long>(xout.size()) ==
        static_cast<amrex::Long>(yout.size()));

    int npts = xout.size();
    for (int i = 0; i < npts; ++i) {
        yout[i] = linear(xinp, yinp, xout[i], ncomp, comp);
    }
}

// Angles in yinp should be defined between [-upper_bound, +upper_bound]
// Use cases would be degrees (upper_bound =360) and radians (upper_bound =2pi)
template <typename C1, typename C2>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    typename std::iterator_traits<C2>::value_type
    linear_angle(
        const C1 xbegin,
        const C1 xend,
        const C2 yinp,
        const typename std::iterator_traits<C1>::value_type& xout,
        const typename std::iterator_traits<C1>::value_type& upper_bound)
{
    using DType1 = typename std::iterator_traits<C1>::value_type;
    const auto idx = bisection_search(xbegin, xend, xout);

    if ((idx.lim == Limits::LOWLIM) || (idx.lim == Limits::UPLIM)) {
        return yinp[idx.idx];
    }
    static constexpr DType1 eps = std::numeric_limits<float>::epsilon();
    const int j = idx.idx;
    const auto denom = (xbegin[j + 1] - xbegin[j]);
    /* https://math.stackexchange.com/questions/2144234/interpolating-between-2-angles
     */
    const auto R = (denom > eps) ? ((xout - xbegin[j]) / denom) : 0.0_rt;
    const auto ub = upper_bound;
    const auto hb = 0.5_rt * upper_bound;
    const auto ohb = 1.5_rt * upper_bound;
    return yinp[j] +
           ((std::fmod((std::fmod(yinp[j + 1] - yinp[j], ub) + ohb), ub) - hb) *
            R);
}

template <typename C1, typename C2>
typename C2::value_type linear_angle(
    const C1& xinp,
    const C2& yinp,
    const typename C1::value_type& xout,
    const typename C1::value_type& upper_bound)
{
    return linear_angle(
        xinp.data(), (xinp.data() + xinp.size()), yinp.data(), xout,
        upper_bound);
}

template <typename C1, typename C2>
void linear_angle(
    const C1& xinp,
    const C2& yinp,
    const C1& xout,
    C2& yout,
    const typename C1::value_type& upper_bound)
{
    AMREX_ASSERT(
        static_cast<amrex::Long>(xinp.size()) ==
        static_cast<amrex::Long>(yinp.size()));
    AMREX_ASSERT(
        static_cast<amrex::Long>(xout.size()) ==
        static_cast<amrex::Long>(yout.size()));

    int npts = xout.size();
    for (int i = 0; i < npts; ++i) {
        yout[i] = linear_angle(xinp, yinp, xout[i], upper_bound);
    }
}

template <typename C1, typename C2>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    typename std::iterator_traits<C2>::value_type
    bilinear_impl(
        const C1 xbegin,
        const C1 ybegin,
        const int ny,
        const C2 zinp,
        const typename std::iterator_traits<C1>::value_type& xout,
        const typename std::iterator_traits<C1>::value_type& yout,
        const Index& xidx,
        const Index& yidx)
{
    using DType1 = typename std::iterator_traits<C1>::value_type;

    if ((xidx.lim == Limits::LOWLIM) || (xidx.lim == Limits::UPLIM) ||
        (yidx.lim == Limits::LOWLIM) || (yidx.lim == Limits::UPLIM)) {
        return zinp[(xidx.idx * ny) + yidx.idx];
    }

    const int i = xidx.idx;
    const int j = yidx.idx;

    static constexpr DType1 eps = std::numeric_limits<float>::epsilon();
    const auto xdenom = (xbegin[i + 1] - xbegin[i]);
    const auto xfacR = (xdenom > eps) ? ((xout - xbegin[i]) / xdenom) : 1.0_rt;
    const auto xfacL = static_cast<DType1>(1.0_rt) - xfacR;
    const auto ydenom = (ybegin[j + 1] - ybegin[j]);
    const auto yfacR = (ydenom > eps) ? ((yout - ybegin[j]) / ydenom) : 1.0_rt;
    const auto yfacL = static_cast<DType1>(1.0_rt) - yfacR;

    const auto z00 = zinp[(i * ny) + j];
    const auto z10 = zinp[((i + 1) * ny) + j];
    const auto z01 = zinp[(i * ny) + (j + 1)];
    const auto z11 = zinp[((i + 1) * ny) + (j + 1)];

    return (z00 * xfacL * yfacL) + (z10 * xfacR * yfacL) +
           (z01 * xfacL * yfacR) + (z11 * xfacR * yfacR);
}

template <typename C1, typename C2>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    typename std::iterator_traits<C2>::value_type
    bilinear(
        const C1 xbegin,
        const C1 xend,
        const C1 ybegin,
        const C1 yend,
        const C2 zinp,
        const typename std::iterator_traits<C1>::value_type& xout,
        const typename std::iterator_traits<C1>::value_type& yout)
{
    const auto xidx = bisection_search(xbegin, xend, xout);
    const auto yidx = bisection_search(ybegin, yend, yout);
    return bilinear_impl(
        xbegin, ybegin, static_cast<int>(yend - ybegin), zinp, xout, yout, xidx,
        yidx);
}

template <typename C1, typename C2>
typename C2::value_type bilinear(
    const C1& xinp,
    const C1& yinp,
    const C2& zinp,
    const typename C1::value_type& xout,
    const typename C1::value_type& yout)
{
    AMREX_ALWAYS_ASSERT((xinp.size() * yinp.size()) == zinp.size());
    return bilinear(
        xinp.data(), (xinp.data() + xinp.size()), yinp.data(),
        (yinp.data() + yinp.size()), zinp.data(), xout, yout);
}

} // namespace amr_wind::interp

#endif /* LINEAR_INTERPOLATION_H */
