#include "amr-wind/utilities/tagging/RefinementCriteria.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind {

RefineCriteriaManager::RefineCriteriaManager(CFDSim& sim) : m_sim(sim) {}

void RefineCriteriaManager::initialize()
{
    BL_PROFILE("amr-wind::RefineCriteriaManager::initialize");
    // Labels for different sampler types
    amrex::Vector<std::string> labels;
    {
        amrex::ParmParse pp("tagging");
        pp.queryarr("labels", labels);
    }

    for (const auto& lbl : labels) {
        const std::string key = "tagging." + lbl;
        amrex::ParmParse pp(key);
        std::string stype;
        pp.get("type", stype);

        auto obj = RefinementCriteria::create(stype, m_sim);
        obj->initialize(key);
        m_refiners.emplace_back(std::move(obj));
    }
}

void RefineCriteriaManager::tag_cells(
    int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
{
    BL_PROFILE("amr-wind::RefineCriteriaManager::tag_cells");
    for (const auto& rc : m_refiners) {
        (*rc)(lev, tags, time, ngrow);
        amrex::Gpu::streamSynchronize();
    }
}

} // namespace amr_wind
