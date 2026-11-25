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

        std::string op_str = "or";
        pp.query("operator", op_str);
        tagging::TaggingOperator op =
            tagging::string_to_operator(amrex::toLower(op_str));

        auto obj = RefinementCriteria::create(stype, m_sim);
        obj->initialize(key);
        m_refiners.push_back({std::move(obj), op});
    }
}

void RefineCriteriaManager::tag_cells(
    int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
{
    BL_PROFILE("kynema-sgf::RefineCriteriaManager::tag_cells");
    if (m_refiners.empty()) {
        return;
    }
    amrex::TagBoxArray current(tags.boxArray(), tags.DistributionMap());
    for (const auto& step : m_refiners) {
        if (!step.criterion->active(lev)) {
            continue;
        }
        current.setVal(amrex::TagBox::CLEAR);
        (*step.criterion)(lev, current, time, ngrow);
        const auto op = step.op;
        for (amrex::MFIter mfi(tags); mfi.isValid(); ++mfi) {
            const auto& previous = tags.array(mfi);
            const auto& predicate = current.const_array(mfi);
            amrex::ParallelFor(
                mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    previous(i, j, k) = tagging::tag_val(
                        previous(i, j, k) == amrex::TagBox::SET,
                        predicate(i, j, k) == amrex::TagBox::SET, op);
                });
        }
        amrex::Gpu::streamSynchronize();
    }
}

} // namespace kynema_sgf
