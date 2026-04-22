#include "turb_utils.H"

#include <set>

#include "AMReX_ParmParse.H"

namespace amr_wind::turbulence::turb_utils {

void inject_turbulence_src_terms(
    const std::string& key, const amrex::Vector<std::string>& terms)
{
    amrex::ParmParse pp(key);
    amrex::Vector<std::string> srcin;
    pp.queryarr("source_terms", srcin);

    // If there are no sources (usually the case) then add sources and return
    if (srcin.empty()) {
        pp.addarr("source_terms", terms);
        return;
    }

    // There are source terms defined by user... so need to combine terms
    std::set<std::string> sterm_set;
    for (const auto& sname : srcin) {
        sterm_set.insert(sname);
    }
    for (const auto& sname : terms) {
        sterm_set.insert(sname);
    }

    // Now convert the combined terms back into a vector for insertion in input
    // dictionary
    srcin.clear();

    // Prefer to use emplace_back here
    for (const auto& sname : sterm_set) {
        srcin.emplace_back(sname);
    }
    pp.addarr("source_terms", srcin);
}

} // namespace amr_wind::turbulence::turb_utils
