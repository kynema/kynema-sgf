#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include "gtest/gtest.h"
#include "AMReX_MultiFab.H"
#include "AMReX_Vector.H"
#include "amr-wind/core/FieldRepo.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests::utils {

AMREX_INLINE amrex::Real
field_min(const amr_wind::Field& field, const int icomp = 0)
{
    amrex::Real min_val = std::numeric_limits<amrex::Real>::max();
    for (int lev = 0; lev < field.repo().num_active_levels(); ++lev) {
        min_val = amrex::min<amrex::Real>(min_val, field(lev).min(icomp));
    }

    return min_val;
}

AMREX_INLINE amrex::Real
field_max(const amr_wind::Field& field, const int icomp = 0)
{
    amrex::Real max_val = -std::numeric_limits<amrex::Real>::max();
    for (int lev = 0; lev < field.repo().num_active_levels(); ++lev) {
        max_val = amrex::max<amrex::Real>(max_val, field(lev).max(icomp));
    }

    return max_val;
}

AMREX_INLINE amrex::Real
field_min(const amr_wind::ScratchField& field, const int icomp = 0)
{
    amrex::Real min_val = std::numeric_limits<amrex::Real>::max();
    for (int lev = 0; lev < field.repo().num_active_levels(); ++lev) {
        min_val = amrex::min<amrex::Real>(min_val, field(lev).min(icomp));
    }

    return min_val;
}

AMREX_INLINE amrex::Real
field_max(const amr_wind::ScratchField& field, const int icomp = 0)
{
    amrex::Real max_val = -std::numeric_limits<amrex::Real>::max();
    for (int lev = 0; lev < field.repo().num_active_levels(); ++lev) {
        max_val = amrex::max<amrex::Real>(max_val, field(lev).max(icomp));
    }

    return max_val;
}

/** Return the minimum value of a field over all AMR levels
 *
 *  \param nlevels Number of levels that existing in this AMR hierarchy
 *  \param field MultiFabs for the levels
 *  \param ncomp The field component for which min value is desired (default: 0)
 */
AMREX_INLINE amrex::Real field_min(
    const int nlevels,
    const amrex::Vector<amrex::MultiFab*>& field,
    const int ncomp = 0)
{
    amrex::Real min_val = std::numeric_limits<amrex::Real>::max();
    for (int lev = 0; lev < nlevels; ++lev) {
        min_val = amrex::min<amrex::Real>(min_val, field[lev]->min(ncomp));
    }
    return min_val;
}

/** Return the maximum value of a field over all AMR levels
 *
 *  \param nlevels Number of levels that existing in this AMR hierarchy
 *  \param field MultiFabs for the levels
 *  \param ncomp The field component for which max value is desired (default: 0)
 */
AMREX_INLINE amrex::Real field_max(
    const int nlevels,
    const amrex::Vector<amrex::MultiFab*>& field,
    const int ncomp = 0)
{
    amrex::Real max_val = -std::numeric_limits<amrex::Real>::max();
    for (int lev = 0; lev < nlevels; ++lev) {
        max_val = amrex::max<amrex::Real>(max_val, field[lev]->max(ncomp));
    }
    return max_val;
}

/** Return the minimum value all components of a field over all AMR levels
 *
 *  \param nlevels Number of levels that existing in this AMR hierarchy
 *  \param field MultiFabs for the levels
 *  \param min_val Array of component values that are populated
 */
AMREX_INLINE void field_min(
    const int nlevels,
    const amrex::Vector<amrex::MultiFab*>& field,
    amrex::Vector<amrex::Real>& min_val)
{
    for (int icomp = 0; icomp < min_val.size(); ++icomp) {
        min_val[icomp] = field_min(nlevels, field, icomp);
    }
}

/** Return the maximum value all components of a field over all AMR levels
 *
 *  \param nlevels Number of levels that existing in this AMR hierarchy
 *  \param field MultiFabs for the levels
 *  \param max_val Array of component values that are populated
 */
AMREX_INLINE void field_max(
    const int nlevels,
    const amrex::Vector<amrex::MultiFab*>& field,
    amrex::Vector<amrex::Real>& max_val)
{
    for (int icomp = 0; icomp < max_val.size(); ++icomp) {
        max_val[icomp] = field_max(nlevels, field, icomp);
    }
}

AMREX_INLINE void field_minmax(
    const int nlevels,
    const amrex::Vector<amrex::MultiFab*>& field,
    amrex::Real& min_val,
    amrex::Real& max_val,
    const int ncomp = 0)
{
    min_val = field_min(nlevels, field, ncomp);
    max_val = field_max(nlevels, field, ncomp);
}

AMREX_INLINE void field_minmax(
    const int nlevels,
    const amrex::Vector<amrex::MultiFab*>& field,
    amrex::Vector<amrex::Real>& min_val,
    amrex::Vector<amrex::Real>& max_val)
{
    field_min(nlevels, field, min_val);
    field_max(nlevels, field, max_val);
}

template <typename T, typename FType>
AMREX_INLINE T field_probe_impl(
    const FType& field,
    const int lev,
    const int i,
    const int j,
    const int k,
    const int icomp = 0)
{
    const amrex::IntVect iv_val = {i, j, k};
    T value = amrex::ReduceSum(
        field(lev), 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx,
            amrex::Array4<T const> const& field_arr) -> T {
            T val = 0.0;

            if (bx.contains(iv_val)) {
                amrex::Loop(bx, [=, &val](int i0, int j0, int k0) {
                    const amrex::IntVect iv = {i0, j0, k0};
                    if (iv == iv_val) {
                        val = field_arr(iv, icomp);
                    }
                });
            }
            return val;
        });
    return value;
}

AMREX_INLINE int field_probe(
    const amr_wind::IntField& field,
    const int lev,
    const int i,
    const int j,
    const int k,
    const int icomp = 0)
{
    int value = field_probe_impl<int>(field, lev, i, j, k, icomp);
    amrex::ParallelDescriptor::ReduceIntSum(value);
    return value;
}

AMREX_INLINE amrex::Real field_probe(
    const amr_wind::Field& field,
    const int lev,
    const int i,
    const int j,
    const int k,
    const int icomp = 0)
{
    auto value = field_probe_impl<amrex::Real>(field, lev, i, j, k, icomp);
    amrex::ParallelDescriptor::ReduceRealSum(value);
    return value;
}

} // namespace amr_wind_tests::utils

#endif /* TEST_UTILS_H */
