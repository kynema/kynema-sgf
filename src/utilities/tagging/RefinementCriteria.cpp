#include "src/utilities/tagging/RefinementCriteria.H"
#include "src/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace kynema_sgf {

RefineCriteriaManager::RefineCriteriaManager(CFDSim& sim) : m_sim(sim) {}

void RefineCriteriaManager::initialize()
{
    BL_PROFILE("kynema-sgf::RefineCriteriaManager::initialize");
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
        
        std::string op{obj->tag_operator()};
        pp.query("operator", op);
        op = amrex::toLower(op);
        obj->tag_operator() = op;
        if (op != "and" && op != "or" && op != "and_not" && op != "or_not") {
            amrex::Abort(
                "Invalid operator for tagging criteria: " + lbl +
                ". Must be one of: and, or, and_not, or_not.");
        }
        m_refiners.emplace_back(std::move(obj));
    }
}

void RefineCriteriaManager::tag_cells(
    int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
{
    BL_PROFILE("kynema-sgf::RefineCriteriaManager::tag_cells");
    for (const auto& rc : m_refiners) {
        (*rc)(lev, tags, time, ngrow);
        amrex::Gpu::streamSynchronize();
    }
}

} // namespace kynema_sgf
