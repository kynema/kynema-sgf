#ifndef TRANSVERSE_GHOST_FILL_H
#define TRANSVERSE_GHOST_FILL_H

#include "AMReX_MultiFab.H"
#include "AMReX_Geometry.H"
#include "AMReX_BCRec.H"
#include "AMReX_Orientation.H"

namespace kynema_sgf {

/** Complete an auxiliary-filled physical face across transverse grid, periodic,
 * and physical boundaries. Higher-priority Dirichlet faces retain their values.
 * The face's normal ghost layers must already be filled over valid transverse
 * cell-centered indices. Components in bcrec start at orig_comp; destination
 * components at dcomp.
 */
void fill_transverse_ghosts(
    amrex::MultiFab& mfab,
    const amrex::Geometry& geom,
    const amrex::Vector<amrex::BCRec>& bcrec,
    amrex::Orientation ori,
    int dcomp,
    int ncomp,
    int orig_comp = 0);

} // namespace kynema_sgf
#endif
