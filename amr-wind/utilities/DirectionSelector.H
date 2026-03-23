#ifndef DIRECTIONSELECTOR_H
#define DIRECTIONSELECTOR_H

#include "AMReX_Box.H"
#include "AMReX_Gpu.H"

/** select index based on direction input */
template <int Index>
struct DirectionSelector
{
    int operator()(int i, int j, int k) const;
};

template <>
struct DirectionSelector<0>
{
    AMREX_GPU_DEVICE int operator()(int i, int /*unused*/, int /*unused*/) const
    {
        return i;
    };
    int odir1 = 1;
    int odir2 = 2;
};

template <>
struct DirectionSelector<1>
{
    AMREX_GPU_DEVICE int operator()(int /*unused*/, int j, int /*unused*/) const
    {
        return j;
    };
    int odir1 = 0;
    int odir2 = 2;
};

template <>
struct DirectionSelector<2>
{
    AMREX_GPU_DEVICE int operator()(int /*unused*/, int /*unused*/, int k) const
    {
        return k;
    };
    int odir1 = 0;
    int odir2 = 1;
};

using XDir = DirectionSelector<0>;
using YDir = DirectionSelector<1>;
using ZDir = DirectionSelector<2>;

// Given a box, return a 2D box perpendicular to the selected axis.
// For example, if we're using ZDir, return a box covering the x-y plane.
// The IntVect is used to set the constant index in the parallel direction.
template <typename IndexSelector>
AMREX_GPU_HOST_DEVICE amrex::Box
perpendicular_box(const amrex::Box& bx, const amrex::IntVect& iv)
{
    amrex::IntVect plane_lo, plane_hi;

    if (std::is_same_v<IndexSelector, XDir>) {
        plane_lo = {iv[0], bx.smallEnd(1), bx.smallEnd(2)};
        plane_hi = {iv[0], bx.bigEnd(1), bx.bigEnd(2)};
    } else if (std::is_same_v<IndexSelector, YDir>) {
        plane_lo = {bx.smallEnd(0), iv[1], bx.smallEnd(2)};
        plane_hi = {bx.bigEnd(0), iv[1], bx.bigEnd(2)};
    } else {
        plane_lo = {bx.smallEnd(0), bx.smallEnd(1), iv[2]};
        plane_hi = {bx.bigEnd(0), bx.bigEnd(1), iv[2]};
    }

    amrex::Box pbx(plane_lo, plane_hi);

    return pbx;
}

// Given a box, return a 1D box parallel to the selected axis.
// For example, if we're using ZDir, return a box covering the z axis.
// The IntVect is used to set the constant indices in the perpendicular
// direction.
template <typename IndexSelector>
AMREX_GPU_DEVICE amrex::Box
parallel_box(const amrex::Box& bx, const amrex::IntVect& iv)
{
    amrex::IntVect line_lo, line_hi;

    if (std::is_same_v<IndexSelector, XDir>) {
        line_lo = {bx.smallEnd(0), iv[1], iv[2]};
        line_hi = {bx.bigEnd(0), iv[1], iv[2]};
    } else if (std::is_same_v<IndexSelector, YDir>) {
        line_lo = {iv[0], bx.smallEnd(1), iv[2]};
        line_hi = {iv[0], bx.bigEnd(1), iv[2]};
    } else {
        line_lo = {iv[0], iv[1], bx.smallEnd(2)};
        line_hi = {iv[0], iv[1], bx.bigEnd(2)};
    }

    amrex::Box lbx(line_lo, line_hi);

    return lbx;
}

namespace amr_wind {

AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
direction_selector(const int i, const int j, const int k, const int dir)
{
    return (dir == 0) ? i : (dir == 1) ? j : k;
}

} /* namespace amr_wind */
#endif /* DirectionSelector_H */
