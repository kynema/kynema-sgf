#include "amr-wind/utilities/tagging/GeometryRefinement.H"
#include "amr-wind/CFDSim.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

GeometryRefinement::GeometryRefinement(const CFDSim& sim)
    : m_sim(sim), m_max_level(m_sim.mesh().maxLevel())
{}

void GeometryRefinement::initialize(const std::string& key)
{
    amrex::Vector<std::string> shapes;
    {
        amrex::ParmParse pp(key);
        pp.getarr("shapes", shapes);
        pp.query("level", m_set_level);
        pp.query("min_level", m_min_level);
        pp.query("max_level", m_max_level);
    }

    for (const auto& geom : shapes) {
        const std::string& args = key + "." + geom;
        amrex::ParmParse pp(args);
        std::string gtype;
        pp.get("type", gtype);
        auto obj = tagging::GeometryType::create(gtype, m_sim, args);
        m_geom_refiners.emplace_back(std::move(obj));
    }
}

void GeometryRefinement::operator()(
    int level, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
    // If the user has requested a particular level then check for it and exit
    // early
    if ((m_set_level > -1) && (level != m_set_level)) {
        return;
    }

    // If the user has specified a range of levels, check and return early
    if ((level < m_min_level) || (level > m_max_level)) {
        return;
    }

    const auto& mesh = m_sim.mesh();
    const auto& geom = mesh.Geom(level);

    // We are always guaranteed to have at least one field
    const auto& field_fab = (*m_sim.repo().fields()[0])(level);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(field_fab); mfi.isValid(); ++mfi) {
        const auto& bx = mfi.tilebox();
        const auto& tag = tags.array(mfi);

        for (const auto& gg : m_geom_refiners) {
            const auto& bxa = utils::realbox_to_box(gg->bound_box(), geom);
            const auto& bxi = bx & bxa;
            if (bxi.isEmpty()) {
                continue;
            }

            (*gg)(bx, geom, tag);
        }
    }
}

} // namespace amr_wind
