#ifndef FIELD_OPS_H
#define FIELD_OPS_H

#include "amr-wind/core/Field.H"
#include "AMReX_MultiFab.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

/**
 *  \defgroup field_ops Field operations
 *
 *  Operators for manipulating fields and computing derived quantities.
 *
 *  This group documents functions that provide various operations defined on
 *  the field. These operations can be grouped into two major categories: simple
 *  linear algebra operations (e.g., saxpy), and discretization operators such
 *  as gradient, laplacian computations.
 *
 *  \ingroup fields
 */

namespace amr_wind::field_ops {

/** Add two fields \f$y = x + y\f$
 *  \ingroup field_ops
 *
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] src Field to be added
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
template <typename T1, typename T2>
void add(
    T1& dst,
    const T2& src,
    int srccomp,
    int dstcomp,
    int numcomp,
    const amrex::IntVect& nghost)
{
    const int nlevels = dst.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::MultiFab::Add(
            dst(lev), src(lev), srccomp, dstcomp, numcomp, nghost);
    }
}

/** Add two fields \f$y = x + y\f$
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] src Field to be added
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
template <typename T1, typename T2>
void add(
    T1& dst, const T2& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    add(dst, src, srccomp, dstcomp, numcomp, amrex::IntVect(nghost));
}

/** Divide two fields \f$y = y / x\f$
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] src Field to be divisor
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] ncomp_src Number of components to divide by
 *  \param [in] ncomp_dst Number of components to be divided
 *  \param [in] nghost Number of ghost cells to be updated
 */
template <typename T1, typename T2>
void divide(
    T1& dst,
    const T2& src,
    const int srccomp,
    const int dstcomp,
    const int ncomp_src,
    const int ncomp_dst,
    const amrex::IntVect& nghost)
{
    const int nlevels = dst.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        if (ncomp_dst == ncomp_src) {
            amrex::MultiFab::Divide(
                dst(lev), src(lev), srccomp, dstcomp, ncomp_dst, nghost);
        } else if (ncomp_src == 1) {
            for (int n = 0; n < ncomp_dst; ++n) {
                amrex::MultiFab::Divide(
                    dst(lev), src(lev), srccomp, dstcomp + n, ncomp_src,
                    nghost);
            }
        } else {
            amrex::Abort(
                "field_ops::divide : number of divisor components (ncomp_src) "
                "must either be 1 or match number of dividend components "
                "(ncomp_dst).");
        }
    }
}

/** Divide two fields \f$y = y / x\f$
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] src Field to be divisor
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] ncomp_src Number of components to divide by
 *  \param [in] ncomp_dst Number of components to be divided
 *  \param [in] nghost Number of ghost cells to be updated
 */
template <typename T1, typename T2>
void divide(
    T1& dst,
    const T2& src,
    const int srccomp,
    const int dstcomp,
    const int ncomp_src,
    const int ncomp_dst,
    const int nghost)
{
    divide(
        dst, src, srccomp, dstcomp, ncomp_src, ncomp_dst,
        amrex::IntVect(nghost));
}

/** Copy source field to destination field
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] src Field to be added
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
template <typename T1, typename T2>
void copy(
    T1& dst,
    const T2& src,
    const int srccomp,
    const int dstcomp,
    const int numcomp,
    const amrex::IntVect& nghost)
{
    const int nlevels = dst.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::MultiFab::Copy(
            dst(lev), src(lev), srccomp, dstcomp, numcomp, nghost);
    }
}

/** Copy source field to destination field
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] src Field to be added
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
template <typename T1, typename T2>
void copy(
    T1& dst, const T2& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    copy(dst, src, srccomp, dstcomp, numcomp, amrex::IntVect(nghost));
}

/** Perform operation \f$y = a x + y\f$
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] a  Coefficient
 *  \param [in] src Field to be added
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
template <typename T1, typename T2>
void saxpy(
    T1& dst,
    const amrex::Real a,
    const T2& src,
    const int srccomp,
    const int dstcomp,
    const int numcomp,
    const amrex::IntVect& nghost)
{
    const int nlevels = dst.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::MultiFab::Saxpy(
            dst(lev), a, src(lev), srccomp, dstcomp, numcomp, nghost);
    }
}

/** Perform operation \f$y = a x + y\f$
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] a  Coefficient
 *  \param [in] src Field to be added
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
template <typename T1, typename T2>
void saxpy(
    T1& dst,
    const amrex::Real a,
    const T2& src,
    const int srccomp,
    const int dstcomp,
    const int numcomp,
    const int nghost)
{
    saxpy(dst, a, src, srccomp, dstcomp, numcomp, amrex::IntVect(nghost));
}

/** Perform operation \f$y = x + a y\f$
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] a  Coefficient
 *  \param [in] src Field to be added
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
template <typename T1, typename T2>
void xpay(
    T1& dst,
    const amrex::Real a,
    const T2& src,
    const int srccomp,
    const int dstcomp,
    const int numcomp,
    const amrex::IntVect& nghost)
{
    const int nlevels = dst.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::MultiFab::Xpay(
            dst(lev), a, src(lev), srccomp, dstcomp, numcomp, nghost);
    }
}

/** Perform operation \f$y = x + a y\f$
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] a  Coefficient
 *  \param [in] src Field to be added
 *  \param [in] srccomp Starting component index of source field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
template <typename T1, typename T2>
void xpay(
    T1& dst,
    const amrex::Real a,
    const T2& src,
    const int srccomp,
    const int dstcomp,
    const int numcomp,
    const int nghost)
{
    xpay(dst, a, src, srccomp, dstcomp, numcomp, amrex::IntVect(nghost));
}

/** Perform operation \f$z = a x + b y\f$
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] a  Coefficient
 *  \param [in] x Field to be added
 *  \param [in] xcomp Starting component index of `x` field
 *  \param [in] b  Coefficient
 *  \param [in] y Field to be added
 *  \param [in] ycomp Starting component index of `y` field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
template <typename T1, typename T2, typename T3>
void lincomb(
    T1& dst,
    const amrex::Real a,
    const T2& x,
    const int xcomp,
    const amrex::Real b,
    const T3& y,
    const int ycomp,
    const int dstcomp,
    const int numcomp,
    const amrex::IntVect& nghost)
{
    const int nlevels = dst.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::MultiFab::LinComb(
            dst(lev), a, x(lev), xcomp, b, y(lev), ycomp, dstcomp, numcomp,
            nghost);
    }
}

/** Perform operation \f$z = a x + b y\f$
 *  \ingroup field_ops
 *
 *  \tparam T1 Field or ScratchField
 *  \tparam T2 Field or ScratchField
 *  \param [out] dst Field that is updated
 *  \param [in] a  Coefficient
 *  \param [in] x Field to be added
 *  \param [in] xcomp Starting component index of `x` field
 *  \param [in] b  Coefficient
 *  \param [in] y Field to be added
 *  \param [in] ycomp Starting component index of `y` field
 *  \param [in] dstcomp Starting component index of destination field
 *  \param [in] numcomp Number of components to be updated
 *  \param [in] nghost Number of ghost cells to be updated
 */
template <typename T1, typename T2, typename T3>
void lincomb(
    T1& dst,
    const amrex::Real a,
    const T2& x,
    const int xcomp,
    const amrex::Real b,
    const T3& y,
    const int ycomp,
    const int dstcomp,
    const int numcomp,
    const int nghost)
{
    lincomb(
        dst, a, x, xcomp, b, y, ycomp, dstcomp, numcomp,
        amrex::IntVect(nghost));
}

/** Set the lower bound for a given scalar field
 * \ingroup field_ops
 *
 * \param [in] field Field that will be clipped
 * \param [in] min_value Minimum value for a given field
 * \param [in] icomp Component to set the minimum bound
 */
template <typename FType>
void lower_bound(FType& field, const amrex::Real min_value, const int icomp = 0)
{
    const auto& repo = field.repo();
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& farrs = field(lev).arrays();
        amrex::ParallelFor(
            field(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                farrs[nbx](i, j, k, icomp) =
                    amrex::max(min_value, farrs[nbx](i, j, k, icomp));
            });
    }
    amrex::Gpu::streamSynchronize();
}

/** Computes the global maximum of a field from all levels
 * \ingroup field_ops
 *
 * \param [in] field Field for global maximum magnitude
 */
template <typename FType>
amrex::Real global_max_magnitude(FType& field)
{
    const auto& repo = field.repo();
    const auto ncomp = field.num_comp();

    amrex::Real maxglobal = 0.0_rt;
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        amrex::Real maxglobal_lev = 0.0_rt;
        maxglobal_lev = amrex::ReduceMax(
            field(lev), 0.,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& b,
                amrex::Array4<amrex::Real const> const& field_arr)
                -> amrex::Real {
                amrex::Real mx = 0.0_rt;
                amrex::Loop(b, [=, &mx](int i, int j, int k) {
                    amrex::Real mag = 0.0_rt;
                    for (int icomp = 0; icomp < ncomp; ++icomp) {
                        mag = mag + (field_arr(i, j, k, icomp) *
                                     field_arr(i, j, k, icomp));
                    }
                    mx = amrex::max(std::sqrt(mag), mx);
                });
                return mx;
            });
        maxglobal = amrex::max(maxglobal, maxglobal_lev);
    }

    amrex::ParallelAllReduce::Max<amrex::Real>(
        maxglobal, amrex::ParallelContext::CommunicatorSub());
    return maxglobal;
}

/** Normalizes a field using its magnitude
 * \ingroup field_ops
 *
 * \param [in] field Field that needs to be normalized
 */
template <typename FType>
void normalize(FType& field)
{
    const amrex::Real eps =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    const auto& repo = field.repo();
    const auto ncomp = field.num_comp();

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& farrs = field(lev).arrays();
        amrex::ParallelFor(
            field(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                auto farr = farrs[nbx];
                // Compute magnitude
                amrex::Real mag = 0.0_rt;
                for (int icomp = 0; icomp < ncomp; ++icomp) {
                    mag = mag + (farr(i, j, k, icomp) * farr(i, j, k, icomp));
                }
                if (mag > eps) {
                    for (int icomp = 0; icomp < ncomp; ++icomp) {
                        farr(i, j, k, icomp) =
                            farr(i, j, k, icomp) / std::sqrt(mag);
                    }
                }
            });
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace amr_wind::field_ops

#endif /* FIELD_OPS_H */
