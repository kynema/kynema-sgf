#ifndef BC_OPS_H
#define BC_OPS_H

/** \file bc_ops.H
 *
 *  Boundary condition utilities
 */

#include "AMReX_BCRec.H"

namespace kynema_sgf::utils {

//! Return a pair of bools,
//! value is true if at least one boundary has an ext_dir specified
std::pair<bool, bool> has_extdir(amrex::BCRec const* bcrec, int ncomp, int dir);
std::pair<bool, bool>
has_extdir_or_ho(amrex::BCRec const* bcrec, int ncomp, int dir);

} // namespace kynema_sgf::utils

#endif /* BC_OPS_H */
