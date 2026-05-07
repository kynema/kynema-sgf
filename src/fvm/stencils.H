#ifndef STENCILS_H
#define STENCILS_H

/** \file stencils.H
 *  \brief Second-order finite volume discretization stencil coefficients
 *
 *  The stencils are defined using 9 coefficients corresponding to `i+1`, `i`
 *  and `i-1` respectively in each of the three directions (see table below).
 *
 *  Coefficients for the first derivative
 *  ```
 *  |   | +1  | 0   | -1  |
 *  |---+-----+-----+-----|
 *  | i | c00 | c01 | c02 |
 *  | j | c10 | c11 | c12 |
 *  | k | c20 | c21 | c22 |
 *  ```
 *
 *  Coefficients for the second derivative
 *  ```
 *  |   | +1  | 0   | -1  |
 *  |---+-----+-----+-----|
 *  | i | s00 | s01 | s02 |
 *  | j | s10 | s11 | s12 |
 *  | k | s20 | s21 | s22 |
 *  ```
 *
 *  Coefficients for the filter
 *  ```
 *  |   | +1  | 0   | -1  |
 *  |---+-----+-----+-----|
 *  | i | f00 | f01 | f02 |
 *  | j | f10 | f11 | f12 |
 *  | k | f20 | f21 | f22 |
 *  ```
 */

#include "AMReX_REAL.H"
#include "AMReX_Box.H"
#include "AMReX_Orientation.H"
#include "AMReX_Geometry.H"

using namespace amrex::literals;

namespace kynema_sgf::fvm::stencil {

namespace impl {

/** Return a layer of cells adjacent to a non-periodic lower boundary
 *
 *  Given a box, this function returns a new box that contains the first layer
 *  of cells immediately adjacent to the lower boundary in the given direction
 *  (`idir`).
 *
 *  \param [in] bx The box whose layers are being considered
 *  \param [in] geom The amrex::Geometry object
 *  \param [in] idir Direction in which the first cell adjacent to the boundary
 * is desired
 *
 *  \return A box with a single layer of cells adjacent to the boundary or empty
 */
AMREX_INLINE amrex::Box
box_lo(const amrex::Box& bx, const amrex::Geometry& geom, const int idir)
{
    if (geom.isPeriodic(idir)) {
        return {};
    }
    if (bx.smallEnd(idir) != geom.Domain().smallEnd(idir)) {
        return {};
    }

    amrex::IntVect low(bx.smallEnd());
    amrex::IntVect hi(bx.bigEnd());
    int sm = low[idir];
    low.setVal(idir, sm);
    hi.setVal(idir, sm);
    return {low, hi};
}

/** Return a layer of cells adjacent to a non-periodic upper boundary
 *
 *  Given a box, this function returns a new box that contains the first layer
 *  of cells immediately adjacent to the lower boundary in the given direction
 *  (`idir`).
 *
 *  \param [in] bx The box whose layers are being considered
 *  \param [in] geom The amrex::Geometry object
 *  \param [in] idir Direction in which the first cell adjacent to the boundary
 * is desired
 *
 *  \return A box with a single layer of cells adjacent to the boundary or empty
 */
AMREX_INLINE amrex::Box
box_hi(const amrex::Box& bx, const amrex::Geometry& geom, const int idir)
{
    if (geom.isPeriodic(idir)) {
        return {};
    }
    if (bx.bigEnd(idir) != geom.Domain().bigEnd(idir)) {
        return {};
    }

    amrex::IntVect low(bx.smallEnd());
    amrex::IntVect hi(bx.bigEnd());
    int sm = hi[idir];
    low.setVal(idir, sm);
    hi.setVal(idir, sm);
    return {low, hi};
}
} // namespace impl

/** Second-order FVM stencils for interior cells
 */
struct StencilInterior
{
    // First derivatives
    static constexpr amrex::Real c00 = 0.5_rt;
    static constexpr amrex::Real c01 = 0.0_rt;
    static constexpr amrex::Real c02 = -0.5_rt;
    static constexpr amrex::Real c10 = c00;
    static constexpr amrex::Real c11 = c01;
    static constexpr amrex::Real c12 = c02;
    static constexpr amrex::Real c20 = c00;
    static constexpr amrex::Real c21 = c01;
    static constexpr amrex::Real c22 = c02;

    // Second derivatives
    static constexpr amrex::Real s00 = 1.0_rt;
    static constexpr amrex::Real s01 = -2.0_rt;
    static constexpr amrex::Real s02 = 1.0_rt;
    static constexpr amrex::Real s10 = s00;
    static constexpr amrex::Real s11 = s01;
    static constexpr amrex::Real s12 = s02;
    static constexpr amrex::Real s20 = s00;
    static constexpr amrex::Real s21 = s01;
    static constexpr amrex::Real s22 = s02;

    // Filter
    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.5_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = f00;
    static constexpr amrex::Real f11 = f01;
    static constexpr amrex::Real f12 = f02;
    static constexpr amrex::Real f20 = f00;
    static constexpr amrex::Real f21 = f01;
    static constexpr amrex::Real f22 = f02;

    static amrex::Box
    box(const amrex::Box& bx, const amrex::Geometry& /*unused*/)
    {
        return bx;
    }
};

struct StencilILO
{
    static constexpr int direction = 0;
    static constexpr int side = amrex::Orientation::Side::low;

    static constexpr amrex::Real c00 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = 1.0_rt;
    static constexpr amrex::Real c02 = -4.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = 0.5_rt;
    static constexpr amrex::Real c11 = 0.0_rt;
    static constexpr amrex::Real c12 = -0.5_rt;
    static constexpr amrex::Real c20 = c10;
    static constexpr amrex::Real c21 = c11;
    static constexpr amrex::Real c22 = c12;

    static constexpr amrex::Real s00 = 4.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 8.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = 1.0_rt;          // ayy
    static constexpr amrex::Real s11 = -2.0_rt;         // byy
    static constexpr amrex::Real s12 = 1.0_rt;          // cyy
    static constexpr amrex::Real s20 = s10;             // azz
    static constexpr amrex::Real s21 = s11;             // bzz
    static constexpr amrex::Real s22 = s12;             // czz

    static constexpr amrex::Real f00 = 0.125_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = 0.25_rt;
    static constexpr amrex::Real f11 = 0.5_rt;
    static constexpr amrex::Real f12 = 0.25_rt;
    static constexpr amrex::Real f20 = f10;
    static constexpr amrex::Real f21 = f11;
    static constexpr amrex::Real f22 = f12;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return impl::box_lo(bx, geom, StencilILO::direction);
    }
};

struct StencilJLO
{
    static constexpr int direction = 1;
    static constexpr int side = amrex::Orientation::Side::low;

    static constexpr amrex::Real c00 = 0.5_rt;
    static constexpr amrex::Real c01 = 0.0_rt;
    static constexpr amrex::Real c02 = -0.5_rt;
    static constexpr amrex::Real c10 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c11 = 1.0_rt;
    static constexpr amrex::Real c12 = -4.0_rt / 3.0_rt;
    static constexpr amrex::Real c20 = c00;
    static constexpr amrex::Real c21 = c01;
    static constexpr amrex::Real c22 = c02;

    static constexpr amrex::Real s00 = 1.0_rt;          // axx
    static constexpr amrex::Real s01 = -2.0_rt;         // bxx
    static constexpr amrex::Real s02 = 1.0_rt;          // cxx
    static constexpr amrex::Real s10 = 4.0_rt / 3.0_rt; // ayy
    static constexpr amrex::Real s11 = -4.0_rt;         // byy
    static constexpr amrex::Real s12 = 8.0_rt / 3.0_rt; // cyy
    static constexpr amrex::Real s20 = s00;             // azz
    static constexpr amrex::Real s21 = s01;             // bzz
    static constexpr amrex::Real s22 = s02;             // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.5_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = 0.125_rt;
    static constexpr amrex::Real f11 = 0.625_rt;
    static constexpr amrex::Real f12 = 0.25_rt;
    static constexpr amrex::Real f20 = f00;
    static constexpr amrex::Real f21 = f01;
    static constexpr amrex::Real f22 = f02;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return impl::box_lo(bx, geom, StencilJLO::direction);
    }
};

struct StencilKLO
{
    static constexpr int direction = 2;
    static constexpr int side = amrex::Orientation::Side::low;

    static constexpr amrex::Real c00 = 0.5_rt;
    static constexpr amrex::Real c01 = 0.0_rt;
    static constexpr amrex::Real c02 = -0.5_rt;
    static constexpr amrex::Real c10 = c00;
    static constexpr amrex::Real c11 = c01;
    static constexpr amrex::Real c12 = c02;
    static constexpr amrex::Real c20 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c21 = 1.0_rt;
    static constexpr amrex::Real c22 = -4.0_rt / 3.0_rt;

    static constexpr amrex::Real s00 = 1.0_rt;          // axx
    static constexpr amrex::Real s01 = -2.0_rt;         // bxx
    static constexpr amrex::Real s02 = 1.0_rt;          // cxx
    static constexpr amrex::Real s10 = s00;             // ayy
    static constexpr amrex::Real s11 = s01;             // byy
    static constexpr amrex::Real s12 = s02;             // cyy
    static constexpr amrex::Real s20 = 4.0_rt / 3.0_rt; // azz
    static constexpr amrex::Real s21 = -4.0_rt;         // bzz
    static constexpr amrex::Real s22 = 8.0_rt / 3.0_rt; // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.5_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = f00;
    static constexpr amrex::Real f11 = f01;
    static constexpr amrex::Real f12 = f02;
    static constexpr amrex::Real f20 = 0.125_rt;
    static constexpr amrex::Real f21 = 0.625_rt;
    static constexpr amrex::Real f22 = 0.25_rt;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return impl::box_lo(bx, geom, StencilKLO::direction);
    }
};

struct StencilIHI
{
    static constexpr int direction = 0;
    static constexpr int side = amrex::Orientation::Side::high;

    static constexpr amrex::Real c00 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = -1.0_rt;
    static constexpr amrex::Real c02 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = 0.5_rt;
    static constexpr amrex::Real c11 = 0.0_rt;
    static constexpr amrex::Real c12 = -0.5_rt;
    static constexpr amrex::Real c20 = c10;
    static constexpr amrex::Real c21 = c11;
    static constexpr amrex::Real c22 = c12;

    static constexpr amrex::Real s00 = 8.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 4.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = 1.0_rt;          // ayy
    static constexpr amrex::Real s11 = -2.0_rt;         // byy
    static constexpr amrex::Real s12 = 1.0_rt;          // cyy
    static constexpr amrex::Real s20 = s10;             // azz
    static constexpr amrex::Real s21 = s11;             // bzz
    static constexpr amrex::Real s22 = s12;             // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.125_rt;
    static constexpr amrex::Real f10 = 0.25_rt;
    static constexpr amrex::Real f11 = 0.5_rt;
    static constexpr amrex::Real f12 = 0.25_rt;
    static constexpr amrex::Real f20 = f10;
    static constexpr amrex::Real f21 = f11;
    static constexpr amrex::Real f22 = f12;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return impl::box_hi(bx, geom, StencilIHI::direction);
    }
};

struct StencilJHI
{
    static constexpr int direction = 1;
    static constexpr int side = amrex::Orientation::Side::high;

    static constexpr amrex::Real c00 = 0.5_rt;
    static constexpr amrex::Real c01 = 0.0_rt;
    static constexpr amrex::Real c02 = -0.5_rt;
    static constexpr amrex::Real c10 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c11 = -1.0_rt;
    static constexpr amrex::Real c12 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c20 = c00;
    static constexpr amrex::Real c21 = c01;
    static constexpr amrex::Real c22 = c02;

    static constexpr amrex::Real s00 = 1.0_rt;          // axx
    static constexpr amrex::Real s01 = -2.0_rt;         // bxx
    static constexpr amrex::Real s02 = 1.0_rt;          // cxx
    static constexpr amrex::Real s10 = 8.0_rt / 3.0_rt; // ayy
    static constexpr amrex::Real s11 = -4.0_rt;         // byy
    static constexpr amrex::Real s12 = 4.0_rt / 3.0_rt; // cyy
    static constexpr amrex::Real s20 = s00;             // azz
    static constexpr amrex::Real s21 = s01;             // bzz
    static constexpr amrex::Real s22 = s02;             // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.5_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = 0.25_rt;
    static constexpr amrex::Real f11 = 0.625_rt;
    static constexpr amrex::Real f12 = 0.125_rt;
    static constexpr amrex::Real f20 = f00;
    static constexpr amrex::Real f21 = f01;
    static constexpr amrex::Real f22 = f02;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return impl::box_hi(bx, geom, StencilJHI::direction);
    }
};

struct StencilKHI
{
    static constexpr int direction = 2;
    static constexpr int side = amrex::Orientation::Side::high;

    static constexpr amrex::Real c00 = 0.5_rt;
    static constexpr amrex::Real c01 = 0.0_rt;
    static constexpr amrex::Real c02 = -0.5_rt;
    static constexpr amrex::Real c10 = c00;
    static constexpr amrex::Real c11 = c01;
    static constexpr amrex::Real c12 = c02;
    static constexpr amrex::Real c20 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c21 = -1.0_rt;
    static constexpr amrex::Real c22 = -1.0_rt / 3.0_rt;

    static constexpr amrex::Real s00 = 1.0_rt;          // axx
    static constexpr amrex::Real s01 = -2.0_rt;         // bxx
    static constexpr amrex::Real s02 = 1.0_rt;          // cxx
    static constexpr amrex::Real s10 = s00;             // ayy
    static constexpr amrex::Real s11 = s01;             // byy
    static constexpr amrex::Real s12 = s02;             // cyy
    static constexpr amrex::Real s20 = 8.0_rt / 3.0_rt; // azz
    static constexpr amrex::Real s21 = -4.0_rt;         // bzz
    static constexpr amrex::Real s22 = 4.0_rt / 3.0_rt; // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.5_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = f00;
    static constexpr amrex::Real f11 = f01;
    static constexpr amrex::Real f12 = f02;
    static constexpr amrex::Real f20 = 0.25_rt;
    static constexpr amrex::Real f21 = 0.625_rt;
    static constexpr amrex::Real f22 = 0.125_rt;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return impl::box_hi(bx, geom, StencilKHI::direction);
    }
};

struct StencilIloJlo
{
    static constexpr int direction = 10;

    static constexpr amrex::Real c00 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = 1.0_rt;
    static constexpr amrex::Real c02 = -4.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = c00;
    static constexpr amrex::Real c11 = c01;
    static constexpr amrex::Real c12 = c02;
    static constexpr amrex::Real c20 = 0.5_rt;
    static constexpr amrex::Real c21 = 0.0_rt;
    static constexpr amrex::Real c22 = -0.5_rt;

    static constexpr amrex::Real s00 = 4.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 8.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = s00;             // ayy
    static constexpr amrex::Real s11 = s01;             // byy
    static constexpr amrex::Real s12 = s02;             // cyy
    static constexpr amrex::Real s20 = 1.0_rt;          // azz
    static constexpr amrex::Real s21 = -2.0_rt;         // bzz
    static constexpr amrex::Real s22 = 1.0_rt;          // czz

    static constexpr amrex::Real f00 = 0.125_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = f00;
    static constexpr amrex::Real f11 = f01;
    static constexpr amrex::Real f12 = f02;
    static constexpr amrex::Real f20 = 0.25_rt;
    static constexpr amrex::Real f21 = 0.5_rt;
    static constexpr amrex::Real f22 = 0.25_rt;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilILO::box(bx, geom) & StencilJLO::box(bx, geom);
    }
};

struct StencilJloKlo
{
    static constexpr int direction = 11;

    static constexpr amrex::Real c00 = 0.5_rt;
    static constexpr amrex::Real c01 = 0.0_rt;
    static constexpr amrex::Real c02 = -0.5_rt;
    static constexpr amrex::Real c10 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c11 = 1.0_rt;
    static constexpr amrex::Real c12 = -4.0_rt / 3.0_rt;
    static constexpr amrex::Real c20 = c10;
    static constexpr amrex::Real c21 = c11;
    static constexpr amrex::Real c22 = c12;

    static constexpr amrex::Real s00 = 1.0_rt;          // axx
    static constexpr amrex::Real s01 = -2.0_rt;         // bxx
    static constexpr amrex::Real s02 = 1.0_rt;          // cxx
    static constexpr amrex::Real s10 = 4.0_rt / 3.0_rt; // ayy
    static constexpr amrex::Real s11 = -4.0_rt;         // byy
    static constexpr amrex::Real s12 = 8.0_rt / 3.0_rt; // cyy
    static constexpr amrex::Real s20 = s10;             // azz
    static constexpr amrex::Real s21 = s11;             // bzz
    static constexpr amrex::Real s22 = s12;             // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.5_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = 0.125_rt;
    static constexpr amrex::Real f11 = 0.625_rt;
    static constexpr amrex::Real f12 = 0.25_rt;
    static constexpr amrex::Real f20 = f10;
    static constexpr amrex::Real f21 = f11;
    static constexpr amrex::Real f22 = f12;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilJLO::box(bx, geom) & StencilKLO::box(bx, geom);
    }
};

struct StencilIloKlo
{
    static constexpr int direction = 12;

    static constexpr amrex::Real c00 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = 1.0_rt;
    static constexpr amrex::Real c02 = -4.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = 0.5_rt;
    static constexpr amrex::Real c11 = 0.0_rt;
    static constexpr amrex::Real c12 = -0.5_rt;
    static constexpr amrex::Real c20 = c00;
    static constexpr amrex::Real c21 = c01;
    static constexpr amrex::Real c22 = c02;

    static constexpr amrex::Real s00 = 4.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 8.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = 1.0_rt;          // ayy
    static constexpr amrex::Real s11 = -2.0_rt;         // byy
    static constexpr amrex::Real s12 = 1.0_rt;          // cyy
    static constexpr amrex::Real s20 = s00;             // azz
    static constexpr amrex::Real s21 = s01;             // bzz
    static constexpr amrex::Real s22 = s02;             // czz

    static constexpr amrex::Real f00 = 0.125_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = 0.25_rt;
    static constexpr amrex::Real f11 = 0.5_rt;
    static constexpr amrex::Real f12 = 0.25_rt;
    static constexpr amrex::Real f20 = f00;
    static constexpr amrex::Real f21 = f01;
    static constexpr amrex::Real f22 = f02;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilILO::box(bx, geom) & StencilKLO::box(bx, geom);
    }
};

struct StencilIhiJhi
{
    static constexpr int direction = 13;

    static constexpr amrex::Real c00 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = -1.0_rt;
    static constexpr amrex::Real c02 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = c00;
    static constexpr amrex::Real c11 = c01;
    static constexpr amrex::Real c12 = c02;
    static constexpr amrex::Real c20 = 0.5_rt;
    static constexpr amrex::Real c21 = 0.0_rt;
    static constexpr amrex::Real c22 = -0.5_rt;

    static constexpr amrex::Real s00 = 8.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 4.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = s00;             // ayy
    static constexpr amrex::Real s11 = s01;             // byy
    static constexpr amrex::Real s12 = s02;             // cyy
    static constexpr amrex::Real s20 = 1.0_rt;          // azz
    static constexpr amrex::Real s21 = -2.0_rt;         // bzz
    static constexpr amrex::Real s22 = 1.0_rt;          // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.125_rt;
    static constexpr amrex::Real f10 = f00;
    static constexpr amrex::Real f11 = f01;
    static constexpr amrex::Real f12 = f02;
    static constexpr amrex::Real f20 = 0.25_rt;
    static constexpr amrex::Real f21 = 0.5_rt;
    static constexpr amrex::Real f22 = 0.25_rt;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilIHI::box(bx, geom) & StencilJHI::box(bx, geom);
    }
};

struct StencilJhiKhi
{
    static constexpr int direction = 14;

    static constexpr amrex::Real c00 = 0.5_rt;
    static constexpr amrex::Real c01 = 0.0_rt;
    static constexpr amrex::Real c02 = -0.5_rt;
    static constexpr amrex::Real c10 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c11 = -1.0_rt;
    static constexpr amrex::Real c12 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c20 = c10;
    static constexpr amrex::Real c21 = c11;
    static constexpr amrex::Real c22 = c12;

    static constexpr amrex::Real s00 = 1.0_rt;          // axx
    static constexpr amrex::Real s01 = -2.0_rt;         // bxx
    static constexpr amrex::Real s02 = 1.0_rt;          // cxx
    static constexpr amrex::Real s10 = 8.0_rt / 3.0_rt; // ayy
    static constexpr amrex::Real s11 = -4.0_rt;         // byy
    static constexpr amrex::Real s12 = 4.0_rt / 3.0_rt; // cyy
    static constexpr amrex::Real s20 = s10;             // azz
    static constexpr amrex::Real s21 = s11;             // bzz
    static constexpr amrex::Real s22 = s12;             // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.5_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = 0.25_rt;
    static constexpr amrex::Real f11 = 0.625_rt;
    static constexpr amrex::Real f12 = 0.125_rt;
    static constexpr amrex::Real f20 = f10;
    static constexpr amrex::Real f21 = f11;
    static constexpr amrex::Real f22 = f12;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilJHI::box(bx, geom) & StencilKHI::box(bx, geom);
    }
};

struct StencilIhiKhi
{
    static constexpr int direction = 15;

    static constexpr amrex::Real c00 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = -1.0_rt;
    static constexpr amrex::Real c02 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = 0.5_rt;
    static constexpr amrex::Real c11 = 0.0_rt;
    static constexpr amrex::Real c12 = -0.5_rt;
    static constexpr amrex::Real c20 = c00;
    static constexpr amrex::Real c21 = c01;
    static constexpr amrex::Real c22 = c02;

    static constexpr amrex::Real s00 = 8.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 4.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = 1.0_rt;          // ayy
    static constexpr amrex::Real s11 = -2.0_rt;         // byy
    static constexpr amrex::Real s12 = 1.0_rt;          // cyy
    static constexpr amrex::Real s20 = s00;             // azz
    static constexpr amrex::Real s21 = s01;             // bzz
    static constexpr amrex::Real s22 = s02;             // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.125_rt;
    static constexpr amrex::Real f10 = 0.25_rt;
    static constexpr amrex::Real f11 = 0.5_rt;
    static constexpr amrex::Real f12 = 0.25_rt;
    static constexpr amrex::Real f20 = f00;
    static constexpr amrex::Real f21 = f01;
    static constexpr amrex::Real f22 = f02;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilIHI::box(bx, geom) & StencilKHI::box(bx, geom);
    }
};

struct StencilIloJhi
{
    static constexpr int direction = 16;

    static constexpr amrex::Real c00 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = 1.0_rt;
    static constexpr amrex::Real c02 = -4.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c11 = -1.0_rt;
    static constexpr amrex::Real c12 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c20 = 0.5_rt;
    static constexpr amrex::Real c21 = 0.0_rt;
    static constexpr amrex::Real c22 = -0.5_rt;

    static constexpr amrex::Real s00 = 4.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 8.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = 8.0_rt / 3.0_rt; // ayy
    static constexpr amrex::Real s11 = -4.0_rt;         // byy
    static constexpr amrex::Real s12 = 4.0_rt / 3.0_rt; // cyy
    static constexpr amrex::Real s20 = 1.0_rt;          // azz
    static constexpr amrex::Real s21 = -2.0_rt;         // bzz
    static constexpr amrex::Real s22 = 1.0_rt;          // czz

    static constexpr amrex::Real f00 = 0.125_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = 0.25_rt;
    static constexpr amrex::Real f11 = 0.625_rt;
    static constexpr amrex::Real f12 = 0.125_rt;
    static constexpr amrex::Real f20 = 0.25_rt;
    static constexpr amrex::Real f21 = 0.5_rt;
    static constexpr amrex::Real f22 = 0.25_rt;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilILO::box(bx, geom) & StencilJHI::box(bx, geom);
    }
};

struct StencilJloKhi
{
    static constexpr int direction = 17;

    static constexpr amrex::Real c00 = 0.5_rt;
    static constexpr amrex::Real c01 = 0.0_rt;
    static constexpr amrex::Real c02 = -0.5_rt;
    static constexpr amrex::Real c10 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c11 = 1.0_rt;
    static constexpr amrex::Real c12 = -4.0_rt / 3.0_rt;
    static constexpr amrex::Real c20 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c21 = -1.0_rt;
    static constexpr amrex::Real c22 = -1.0_rt / 3.0_rt;

    static constexpr amrex::Real s00 = 1.0_rt;          // axx
    static constexpr amrex::Real s01 = -2.0_rt;         // bxx
    static constexpr amrex::Real s02 = 1.0_rt;          // cxx
    static constexpr amrex::Real s10 = 4.0_rt / 3.0_rt; // ayy
    static constexpr amrex::Real s11 = -4.0_rt;         // byy
    static constexpr amrex::Real s12 = 8.0_rt / 3.0_rt; // cyy
    static constexpr amrex::Real s20 = 8.0_rt / 3.0_rt; // azz
    static constexpr amrex::Real s21 = -4.0_rt;         // bzz
    static constexpr amrex::Real s22 = 4.0_rt / 3.0_rt; // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.5_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = 0.125_rt;
    static constexpr amrex::Real f11 = 0.625_rt;
    static constexpr amrex::Real f12 = 0.25_rt;
    static constexpr amrex::Real f20 = 0.25_rt;
    static constexpr amrex::Real f21 = 0.625_rt;
    static constexpr amrex::Real f22 = 0.125_rt;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilJLO::box(bx, geom) & StencilKHI::box(bx, geom);
    }
};

struct StencilIloKhi
{
    static constexpr int direction = 18;

    static constexpr amrex::Real c00 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = 1.0_rt;
    static constexpr amrex::Real c02 = -4.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = 0.5_rt;
    static constexpr amrex::Real c11 = 0.0_rt;
    static constexpr amrex::Real c12 = -0.5_rt;
    static constexpr amrex::Real c20 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c21 = -1.0_rt;
    static constexpr amrex::Real c22 = -1.0_rt / 3.0_rt;

    static constexpr amrex::Real s00 = 4.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 8.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = 1;               // ayy
    static constexpr amrex::Real s11 = -2;              // byy
    static constexpr amrex::Real s12 = 1;               // cyy
    static constexpr amrex::Real s20 = 8.0_rt / 3.0_rt; // azz
    static constexpr amrex::Real s21 = -4.0_rt;         // bzz
    static constexpr amrex::Real s22 = 4.0_rt / 3.0_rt; // czz

    static constexpr amrex::Real f00 = 0.125_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = 0.25_rt;
    static constexpr amrex::Real f11 = 0.5_rt;
    static constexpr amrex::Real f12 = 0.25_rt;
    static constexpr amrex::Real f20 = 0.25_rt;
    static constexpr amrex::Real f21 = 0.625_rt;
    static constexpr amrex::Real f22 = 0.125_rt;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilILO::box(bx, geom) & StencilKHI::box(bx, geom);
    }
};

struct StencilIhiJlo
{
    static constexpr int direction = 19;

    static constexpr amrex::Real c00 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = -1.0_rt;
    static constexpr amrex::Real c02 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c11 = 1;
    static constexpr amrex::Real c12 = -4.0_rt / 3.0_rt;
    static constexpr amrex::Real c20 = 0.5_rt;
    static constexpr amrex::Real c21 = 0.0_rt;
    static constexpr amrex::Real c22 = -0.5_rt;

    static constexpr amrex::Real s00 = 8.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 4.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = 4.0_rt / 3.0_rt; // ayy
    static constexpr amrex::Real s11 = -4.0_rt;         // byy
    static constexpr amrex::Real s12 = 8.0_rt / 3.0_rt; // cyy
    static constexpr amrex::Real s20 = 1.0_rt;          // azz
    static constexpr amrex::Real s21 = -2.0_rt;         // bzz
    static constexpr amrex::Real s22 = 1.0_rt;          // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.125_rt;
    static constexpr amrex::Real f10 = 0.125_rt;
    static constexpr amrex::Real f11 = 0.625_rt;
    static constexpr amrex::Real f12 = 0.25_rt;
    static constexpr amrex::Real f20 = 0.25_rt;
    static constexpr amrex::Real f21 = 0.5_rt;
    static constexpr amrex::Real f22 = 0.25_rt;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilIHI::box(bx, geom) & StencilJLO::box(bx, geom);
    }
};

struct StencilJhiKlo
{
    static constexpr int direction = 20;

    static constexpr amrex::Real c00 = 0.5_rt;
    static constexpr amrex::Real c01 = 0.0_rt;
    static constexpr amrex::Real c02 = -0.5_rt;
    static constexpr amrex::Real c10 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c11 = -1.0_rt;
    static constexpr amrex::Real c12 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c20 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c21 = 1.0_rt;
    static constexpr amrex::Real c22 = -4.0_rt / 3.0_rt;

    static constexpr amrex::Real s00 = 1.0_rt;          // axx
    static constexpr amrex::Real s01 = -2.0_rt;         // bxx
    static constexpr amrex::Real s02 = 1.0_rt;          // cxx
    static constexpr amrex::Real s10 = 8.0_rt / 3.0_rt; // ayy
    static constexpr amrex::Real s11 = -4.0_rt;         // byy
    static constexpr amrex::Real s12 = 4.0_rt / 3.0_rt; // cyy
    static constexpr amrex::Real s20 = 4.0_rt / 3.0_rt; // azz
    static constexpr amrex::Real s21 = -4.0_rt;         // bzz
    static constexpr amrex::Real s22 = 8.0_rt / 3.0_rt; // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.5_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = 0.25_rt;
    static constexpr amrex::Real f11 = 0.625_rt;
    static constexpr amrex::Real f12 = 0.125_rt;
    static constexpr amrex::Real f20 = 0.125_rt;
    static constexpr amrex::Real f21 = 0.625_rt;
    static constexpr amrex::Real f22 = 0.25_rt;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilJHI::box(bx, geom) & StencilKLO::box(bx, geom);
    }
};

struct StencilIhiKlo
{
    static constexpr int direction = 21;

    static constexpr amrex::Real c00 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = -1.0_rt;
    static constexpr amrex::Real c02 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = 0.5_rt;
    static constexpr amrex::Real c11 = 0.0_rt;
    static constexpr amrex::Real c12 = -0.5_rt;
    static constexpr amrex::Real c20 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c21 = 1.0_rt;
    static constexpr amrex::Real c22 = -4.0_rt / 3.0_rt;

    static constexpr amrex::Real s00 = 8.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 4.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = 1.0_rt;          // ayy
    static constexpr amrex::Real s11 = -2.0_rt;         // byy
    static constexpr amrex::Real s12 = 1.0_rt;          // cyy
    static constexpr amrex::Real s20 = 4.0_rt / 3.0_rt; // azz
    static constexpr amrex::Real s21 = -4.0_rt;         // bzz
    static constexpr amrex::Real s22 = 8.0_rt / 3.0_rt; // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.125_rt;
    static constexpr amrex::Real f10 = 0.25_rt;
    static constexpr amrex::Real f11 = 0.5_rt;
    static constexpr amrex::Real f12 = 0.25_rt;
    static constexpr amrex::Real f20 = 0.125_rt;
    static constexpr amrex::Real f21 = 0.625_rt;
    static constexpr amrex::Real f22 = 0.25_rt;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilIHI::box(bx, geom) & StencilKLO::box(bx, geom);
    }
};

struct StencilIloJloKlo
{
    static constexpr int direction = 22;

    static constexpr amrex::Real c00 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = 1.0_rt;
    static constexpr amrex::Real c02 = -4.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = c00;
    static constexpr amrex::Real c11 = c01;
    static constexpr amrex::Real c12 = c02;
    static constexpr amrex::Real c20 = c00;
    static constexpr amrex::Real c21 = c01;
    static constexpr amrex::Real c22 = c02;

    static constexpr amrex::Real s00 = 4.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 8.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = s00;             // ayy
    static constexpr amrex::Real s11 = s01;             // byy
    static constexpr amrex::Real s12 = s02;             // cyy
    static constexpr amrex::Real s20 = s00;             // azz
    static constexpr amrex::Real s21 = s01;             // bzz
    static constexpr amrex::Real s22 = s02;             // czz

    static constexpr amrex::Real f00 = 0.125_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = f00;
    static constexpr amrex::Real f11 = f01;
    static constexpr amrex::Real f12 = f02;
    static constexpr amrex::Real f20 = f00;
    static constexpr amrex::Real f21 = f01;
    static constexpr amrex::Real f22 = f02;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilILO::box(bx, geom) & StencilJLO::box(bx, geom) &
               StencilKLO::box(bx, geom);
    }
};

struct StencilIloJhiKlo
{
    static constexpr int direction = 23;

    static constexpr amrex::Real c00 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = 1.0_rt;
    static constexpr amrex::Real c02 = -4.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c11 = -1.0_rt;
    static constexpr amrex::Real c12 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c20 = c00;
    static constexpr amrex::Real c21 = c01;
    static constexpr amrex::Real c22 = c02;

    static constexpr amrex::Real s00 = 4.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 8.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = 8.0_rt / 3.0_rt; // ayy
    static constexpr amrex::Real s11 = -4.0_rt;         // byy
    static constexpr amrex::Real s12 = 4.0_rt / 3.0_rt; // cyy
    static constexpr amrex::Real s20 = s00;             // azz
    static constexpr amrex::Real s21 = s01;             // bzz
    static constexpr amrex::Real s22 = s02;             // czz

    static constexpr amrex::Real f00 = 0.125_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = 0.25_rt;
    static constexpr amrex::Real f11 = 0.625_rt;
    static constexpr amrex::Real f12 = 0.125_rt;
    static constexpr amrex::Real f20 = f00;
    static constexpr amrex::Real f21 = f01;
    static constexpr amrex::Real f22 = f02;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilILO::box(bx, geom) & StencilJHI::box(bx, geom) &
               StencilKLO::box(bx, geom);
    }
};

struct StencilIhiJloKlo
{
    static constexpr int direction = 24;

    static constexpr amrex::Real c00 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = -1.0_rt;
    static constexpr amrex::Real c02 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c11 = 1.0_rt;
    static constexpr amrex::Real c12 = -4.0_rt / 3.0_rt;
    static constexpr amrex::Real c20 = c10;
    static constexpr amrex::Real c21 = c11;
    static constexpr amrex::Real c22 = c12;

    static constexpr amrex::Real s00 = 8.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 4.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = 4.0_rt / 3.0_rt; // ayy
    static constexpr amrex::Real s11 = -4.0_rt;         // byy
    static constexpr amrex::Real s12 = 8.0_rt / 3.0_rt; // cyy
    static constexpr amrex::Real s20 = s10;             // azz
    static constexpr amrex::Real s21 = s11;             // bzz
    static constexpr amrex::Real s22 = s12;             // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.125_rt;
    static constexpr amrex::Real f10 = 0.125_rt;
    static constexpr amrex::Real f11 = 0.625_rt;
    static constexpr amrex::Real f12 = 0.25_rt;
    static constexpr amrex::Real f20 = f10;
    static constexpr amrex::Real f21 = f11;
    static constexpr amrex::Real f22 = f12;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilIHI::box(bx, geom) & StencilJLO::box(bx, geom) &
               StencilKLO::box(bx, geom);
    }
};

struct StencilIhiJhiKlo
{
    static constexpr int direction = 25;

    static constexpr amrex::Real c00 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = -1.0_rt;
    static constexpr amrex::Real c02 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c11 = -1.0_rt;
    static constexpr amrex::Real c12 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c20 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c21 = 1.0_rt;
    static constexpr amrex::Real c22 = -4.0_rt / 3.0_rt;

    static constexpr amrex::Real s00 = 8.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 4.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = 8.0_rt / 3.0_rt; // ayy
    static constexpr amrex::Real s11 = -4.0_rt;         // byy
    static constexpr amrex::Real s12 = 4.0_rt / 3.0_rt; // cyy
    static constexpr amrex::Real s20 = 4.0_rt / 3.0_rt; // azz
    static constexpr amrex::Real s21 = -4.0_rt;         // bzz
    static constexpr amrex::Real s22 = 8.0_rt / 3.0_rt; // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.125_rt;
    static constexpr amrex::Real f10 = 0.25_rt;
    static constexpr amrex::Real f11 = 0.625_rt;
    static constexpr amrex::Real f12 = 0.125_rt;
    static constexpr amrex::Real f20 = 0.125_rt;
    static constexpr amrex::Real f21 = 0.625_rt;
    static constexpr amrex::Real f22 = 0.25_rt;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilIHI::box(bx, geom) & StencilJHI::box(bx, geom) &
               StencilKLO::box(bx, geom);
    }
};

struct StencilIloJloKhi
{
    static constexpr int direction = 26;

    static constexpr amrex::Real c00 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = 1.0_rt;
    static constexpr amrex::Real c02 = -4.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = c00;
    static constexpr amrex::Real c11 = c01;
    static constexpr amrex::Real c12 = c02;
    static constexpr amrex::Real c20 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c21 = -1.0_rt;
    static constexpr amrex::Real c22 = -1.0_rt / 3.0_rt;

    static constexpr amrex::Real s00 = 4.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 8.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = s00;             // ayy
    static constexpr amrex::Real s11 = s01;             // byy
    static constexpr amrex::Real s12 = s02;             // cyy
    static constexpr amrex::Real s20 = 8.0_rt / 3.0_rt; // azz
    static constexpr amrex::Real s21 = -4.0_rt;         // bzz
    static constexpr amrex::Real s22 = 4.0_rt / 3.0_rt; // czz

    static constexpr amrex::Real f00 = 0.125_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = f00;
    static constexpr amrex::Real f11 = f01;
    static constexpr amrex::Real f12 = f02;
    static constexpr amrex::Real f20 = 0.25_rt;
    static constexpr amrex::Real f21 = 0.625_rt;
    static constexpr amrex::Real f22 = 0.125_rt;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilILO::box(bx, geom) & StencilJLO::box(bx, geom) &
               StencilKHI::box(bx, geom);
    }
};

struct StencilIloJhiKhi
{
    static constexpr int direction = 27;

    static constexpr amrex::Real c00 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = 1.0_rt;
    static constexpr amrex::Real c02 = -4.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c11 = -1.0_rt;
    static constexpr amrex::Real c12 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c20 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c21 = -1.0_rt;
    static constexpr amrex::Real c22 = -1.0_rt / 3.0_rt;

    static constexpr amrex::Real s00 = 4.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 8.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = 8.0_rt / 3.0_rt; // ayy
    static constexpr amrex::Real s11 = -4.0_rt;         // byy
    static constexpr amrex::Real s12 = 4.0_rt / 3.0_rt; // cyy
    static constexpr amrex::Real s20 = 8.0_rt / 3.0_rt; // azz
    static constexpr amrex::Real s21 = -4.0_rt;         // bzz
    static constexpr amrex::Real s22 = 4.0_rt / 3.0_rt; // czz

    static constexpr amrex::Real f00 = 0.125_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.25_rt;
    static constexpr amrex::Real f10 = 0.25_rt;
    static constexpr amrex::Real f11 = 0.625_rt;
    static constexpr amrex::Real f12 = 0.125_rt;
    static constexpr amrex::Real f20 = 0.25_rt;
    static constexpr amrex::Real f21 = 0.625_rt;
    static constexpr amrex::Real f22 = 0.125_rt;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilILO::box(bx, geom) & StencilJHI::box(bx, geom) &
               StencilKHI::box(bx, geom);
    }
};

struct StencilIhiJloKhi
{
    static constexpr int direction = 28;

    static constexpr amrex::Real c00 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = -1.0_rt;
    static constexpr amrex::Real c02 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = 1.0_rt / 3.0_rt;
    static constexpr amrex::Real c11 = 1.0_rt;
    static constexpr amrex::Real c12 = -4.0_rt / 3.0_rt;
    static constexpr amrex::Real c20 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c21 = -1.0_rt;
    static constexpr amrex::Real c22 = -1.0_rt / 3.0_rt;

    static constexpr amrex::Real s00 = 8.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 4.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = 4.0_rt / 3.0_rt; // ayy
    static constexpr amrex::Real s11 = -4.0_rt;         // byy
    static constexpr amrex::Real s12 = 8.0_rt / 3.0_rt; // cyy
    static constexpr amrex::Real s20 = 8.0_rt / 3.0_rt; // azz
    static constexpr amrex::Real s21 = -4.0_rt;         // bzz
    static constexpr amrex::Real s22 = 4.0_rt / 3.0_rt; // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.125_rt;
    static constexpr amrex::Real f10 = 0.125_rt;
    static constexpr amrex::Real f11 = 0.625_rt;
    static constexpr amrex::Real f12 = 0.25_rt;
    static constexpr amrex::Real f20 = 0.25_rt;
    static constexpr amrex::Real f21 = 0.625_rt;
    static constexpr amrex::Real f22 = 0.125_rt;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilIHI::box(bx, geom) & StencilJLO::box(bx, geom) &
               StencilKHI::box(bx, geom);
    }
};

struct StencilIhiJhiKhi
{
    static constexpr int direction = 29;

    static constexpr amrex::Real c00 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c01 = -1.0_rt;
    static constexpr amrex::Real c02 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c10 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c11 = -1.0_rt;
    static constexpr amrex::Real c12 = -1.0_rt / 3.0_rt;
    static constexpr amrex::Real c20 = 4.0_rt / 3.0_rt;
    static constexpr amrex::Real c21 = -1.0_rt;
    static constexpr amrex::Real c22 = -1.0_rt / 3.0_rt;

    static constexpr amrex::Real s00 = 8.0_rt / 3.0_rt; // axx
    static constexpr amrex::Real s01 = -4.0_rt;         // bxx
    static constexpr amrex::Real s02 = 4.0_rt / 3.0_rt; // cxx
    static constexpr amrex::Real s10 = 8.0_rt / 3.0_rt; // ayy
    static constexpr amrex::Real s11 = -4.0_rt;         // byy
    static constexpr amrex::Real s12 = 4.0_rt / 3.0_rt; // cyy
    static constexpr amrex::Real s20 = 8.0_rt / 3.0_rt; // azz
    static constexpr amrex::Real s21 = -4.0_rt;         // bzz
    static constexpr amrex::Real s22 = 4.0_rt / 3.0_rt; // czz

    static constexpr amrex::Real f00 = 0.25_rt;
    static constexpr amrex::Real f01 = 0.625_rt;
    static constexpr amrex::Real f02 = 0.125_rt;
    static constexpr amrex::Real f10 = 0.25_rt;
    static constexpr amrex::Real f11 = 0.625_rt;
    static constexpr amrex::Real f12 = 0.125_rt;
    static constexpr amrex::Real f20 = 0.25_rt;
    static constexpr amrex::Real f21 = 0.625_rt;
    static constexpr amrex::Real f22 = 0.125_rt;

    static amrex::Box box(const amrex::Box& bx, const amrex::Geometry& geom)
    {
        return StencilIHI::box(bx, geom) & StencilJHI::box(bx, geom) &
               StencilKHI::box(bx, geom);
    }
};

/** Runtime stencil coefficients for use in fused single-kernel dispatch.
 *
 *  All 27 compile-time stencil types (interior, 6 faces, 12 edges, 8 corners)
 *  can be expressed as a combination of three independent per-direction
 *  coefficient sets.  This struct holds those coefficients as runtime values so
 *  a single kernel can select the correct stencil per-cell.
 *
 *  Coefficient layout mirrors the static stencil structs:
 *  - c0x / s0x / f0x : direction 0 (i)
 *  - c1x / s1x / f1x : direction 1 (j)
 *  - c2x / s2x / f2x : direction 2 (k)
 */
struct StencilCoefficients
{
    // First-derivative (gradient / divergence / vorticity / …) coefficients
    amrex::Real c00, c01, c02; ///< direction 0 (i): phi(i+1), phi(i), phi(i-1)
    amrex::Real c10, c11, c12; ///< direction 1 (j)
    amrex::Real c20, c21, c22; ///< direction 2 (k)

    // Second-derivative (laplacian / curvature) coefficients
    amrex::Real s00, s01, s02; ///< direction 0 (i)
    amrex::Real s10, s11, s12; ///< direction 1 (j)
    amrex::Real s20, s21, s22; ///< direction 2 (k)

    // Filter coefficients
    amrex::Real f00, f01, f02; ///< direction 0 (i)
    amrex::Real f10, f11, f12; ///< direction 1 (j)
    amrex::Real f20, f21, f22; ///< direction 2 (k)
};

/** Select stencil coefficients for cell (i, j, k) at runtime.
 *
 *  The stencil is fully separable: the coefficients for each direction depend
 *  only on whether that direction touches a non-periodic domain boundary.
 *
 *  \param [in] i, j, k    Cell indices
 *  \param [in] dlo        Domain lower-corner (inclusive)
 *  \param [in] dhi        Domain upper-corner (inclusive)
 *  \param [in] periodic   Per-direction periodicity flags (non-zero = periodic)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE StencilCoefficients pick_stencil(
    int i,
    int j,
    int k,
    amrex::IntVect const& dlo,
    amrex::IntVect const& dhi,
    amrex::GpuArray<int, AMREX_SPACEDIM> const& periodic) noexcept
{
    StencilCoefficients sc;

    // --- direction 0 (i) ---
    if (!periodic[0] && i == dlo[0]) {
        // ILO: forward-biased one-sided difference
        sc.c00 = 1.0_rt / 3.0_rt;
        sc.c01 = 1.0_rt;
        sc.c02 = -4.0_rt / 3.0_rt;
        sc.s00 = 4.0_rt / 3.0_rt;
        sc.s01 = -4.0_rt;
        sc.s02 = 8.0_rt / 3.0_rt;
        sc.f00 = 0.125_rt;
        sc.f01 = 0.625_rt;
        sc.f02 = 0.25_rt;
    } else if (!periodic[0] && i == dhi[0]) {
        // IHI: backward-biased one-sided difference
        sc.c00 = 4.0_rt / 3.0_rt;
        sc.c01 = -1.0_rt;
        sc.c02 = -1.0_rt / 3.0_rt;
        sc.s00 = 8.0_rt / 3.0_rt;
        sc.s01 = -4.0_rt;
        sc.s02 = 4.0_rt / 3.0_rt;
        sc.f00 = 0.25_rt;
        sc.f01 = 0.625_rt;
        sc.f02 = 0.125_rt;
    } else {
        // Interior: central difference
        sc.c00 = 0.5_rt;
        sc.c01 = 0.0_rt;
        sc.c02 = -0.5_rt;
        sc.s00 = 1.0_rt;
        sc.s01 = -2.0_rt;
        sc.s02 = 1.0_rt;
        sc.f00 = 0.25_rt;
        sc.f01 = 0.5_rt;
        sc.f02 = 0.25_rt;
    }

    // --- direction 1 (j) ---
    if (!periodic[1] && j == dlo[1]) {
        sc.c10 = 1.0_rt / 3.0_rt;
        sc.c11 = 1.0_rt;
        sc.c12 = -4.0_rt / 3.0_rt;
        sc.s10 = 4.0_rt / 3.0_rt;
        sc.s11 = -4.0_rt;
        sc.s12 = 8.0_rt / 3.0_rt;
        sc.f10 = 0.125_rt;
        sc.f11 = 0.625_rt;
        sc.f12 = 0.25_rt;
    } else if (!periodic[1] && j == dhi[1]) {
        sc.c10 = 4.0_rt / 3.0_rt;
        sc.c11 = -1.0_rt;
        sc.c12 = -1.0_rt / 3.0_rt;
        sc.s10 = 8.0_rt / 3.0_rt;
        sc.s11 = -4.0_rt;
        sc.s12 = 4.0_rt / 3.0_rt;
        sc.f10 = 0.25_rt;
        sc.f11 = 0.625_rt;
        sc.f12 = 0.125_rt;
    } else {
        sc.c10 = 0.5_rt;
        sc.c11 = 0.0_rt;
        sc.c12 = -0.5_rt;
        sc.s10 = 1.0_rt;
        sc.s11 = -2.0_rt;
        sc.s12 = 1.0_rt;
        sc.f10 = 0.25_rt;
        sc.f11 = 0.5_rt;
        sc.f12 = 0.25_rt;
    }

    // --- direction 2 (k) ---
    if (!periodic[2] && k == dlo[2]) {
        sc.c20 = 1.0_rt / 3.0_rt;
        sc.c21 = 1.0_rt;
        sc.c22 = -4.0_rt / 3.0_rt;
        sc.s20 = 4.0_rt / 3.0_rt;
        sc.s21 = -4.0_rt;
        sc.s22 = 8.0_rt / 3.0_rt;
        sc.f20 = 0.125_rt;
        sc.f21 = 0.625_rt;
        sc.f22 = 0.25_rt;
    } else if (!periodic[2] && k == dhi[2]) {
        sc.c20 = 4.0_rt / 3.0_rt;
        sc.c21 = -1.0_rt;
        sc.c22 = -1.0_rt / 3.0_rt;
        sc.s20 = 8.0_rt / 3.0_rt;
        sc.s21 = -4.0_rt;
        sc.s22 = 4.0_rt / 3.0_rt;
        sc.f20 = 0.25_rt;
        sc.f21 = 0.625_rt;
        sc.f22 = 0.125_rt;
    } else {
        sc.c20 = 0.5_rt;
        sc.c21 = 0.0_rt;
        sc.c22 = -0.5_rt;
        sc.s20 = 1.0_rt;
        sc.s21 = -2.0_rt;
        sc.s22 = 1.0_rt;
        sc.f20 = 0.25_rt;
        sc.f21 = 0.5_rt;
        sc.f22 = 0.25_rt;
    }

    return sc;
}

} // namespace kynema_sgf::fvm::stencil

#endif /* STENCILS_H */
