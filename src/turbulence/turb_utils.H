#ifndef TURB_UTILS_H
#define TURB_UTILS_H

#include <string>
#include "AMReX_Vector.H"

namespace amr_wind::turbulence::turb_utils {

/** Utility function to inject extra source terms into relevant equations
 *
 *  \param key Name of the PDE
 *  \param terms List of identifiers for the source terms to be added
 */
void inject_turbulence_src_terms(
    const std::string& key, const amrex::Vector<std::string>& terms);

} // namespace amr_wind::turbulence::turb_utils

#endif /* TURB_UTILS_H */
