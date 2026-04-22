#ifndef ITER_TOOLS_H
#define ITER_TOOLS_H

#include "AMReX_MultiFab.H"
#include "AMReX_Vector.H"
#include "src/core/FieldRepo.H"

/** \file iter_tools.H
 *  Iteration utilities
 */

namespace kynema_sgf_tests {

/** Helper to run a functor over all levels of a given field
 *
 *  The functor is a void function that accepts two arguments: the current
 *  level, and an MFIter instance for the MultiFab at the given level.
 *
 *  \param nlevels Number of valid levels in the AMR hierarchy
 *  \param field MultiFabs for each of the `nlevels` levels
 *  \param func Functor to execute over the MultiFabs
 */

template <typename Functor>
void run_algorithm_single_lev(
    const int lev, amrex::MultiFab& lfab, const Functor& func)
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(lfab); mfi.isValid(); ++mfi) {
        func(lev, mfi);
    }
}

template <typename Functor>
void run_algorithm_single_lev(
    const int lev, amrex::iMultiFab& lfab, const Functor& func)
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(lfab); mfi.isValid(); ++mfi) {
        func(lev, mfi);
    }
}

template <typename Functor>
void run_algorithm(
    const int nlevels,
    amrex::Vector<amrex::MultiFab*>& field,
    const Functor& func)
{
    for (int lev = 0; lev < nlevels; ++lev) {
        auto& lfab = *field[lev];
        run_algorithm_single_lev(lev, lfab, func);
    }
}

template <typename Functor>
void run_algorithm(kynema_sgf::Field& field, const Functor& func)
{
    const int nlevels = field.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        auto& lfab = field(lev);
        run_algorithm_single_lev(lev, lfab, func);
    }
}

template <typename Functor>
void run_algorithm(
    const int nlevels,
    amrex::Vector<amrex::iMultiFab*>& intfield,
    const Functor& func)
{
    for (int lev = 0; lev < nlevels; ++lev) {
        auto& lfab = *intfield[lev];
        run_algorithm_single_lev(lev, lfab, func);
    }
}

template <typename Functor>
void run_algorithm(kynema_sgf::IntField& intfield, const Functor& func)
{
    const int nlevels = intfield.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        auto& lfab = intfield(lev);
        run_algorithm_single_lev(lev, lfab, func);
    }
}

} // namespace kynema_sgf_tests

#endif /* ITER_TOOLS_H */
