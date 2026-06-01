#ifndef PP_UTILS_H
#define PP_UTILS_H

#include <string>
#include "AMReX_ParmParse.H"
#include "AMReX_Vector.H"

namespace kynema_sgf_tests::pp_utils {

template <typename T>
void set_default(amrex::ParmParse& pp, const std::string& key, const T& val)
{
    if (!pp.contains(key)) {
        pp.add(key.data(), val);
    }
}

template <typename T>
void set_default(
    amrex::ParmParse& pp, const std::string& key, const amrex::Vector<T>& val)
{
    if (!pp.contains(key)) {
        pp.addarr(key.data(), val);
    }
}

//! Default inputs for kynema_sgf::SimTime
void default_time_inputs();

//! Default inputs for amrex::AmrCore
void default_mesh_inputs();

//! Flag indicating whether managed memory is enabled
bool has_managed_memory();

} // namespace kynema_sgf_tests::pp_utils

#endif /* PP_UTILS_H */
