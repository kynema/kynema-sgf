#include "amr-wind/utilities/tagging/OversetRefinement.H"
#include "amr-wind/CFDSim.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {

OversetRefinement::OversetRefinement(const CFDSim& sim)
    : m_sim(sim), m_max_lev(m_sim.mesh().maxLevel())
{}

void OversetRefinement::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);
    pp.query("max_level", m_max_lev);
    pp.query("tag_fringe", m_tag_fringe);
    pp.query("tag_hole", m_tag_hole);
}

void OversetRefinement::operator()(
    int level, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
    if (level > m_max_lev) {
        return;
    }

    const auto& ibcell = m_sim.repo().get_int_field("iblank_cell");
    const auto& ibfab = ibcell(level);

    const bool tag_fringe = m_tag_fringe;
    const bool tag_hole = m_tag_hole;

    const auto& ibarrs = ibfab.const_arrays();
    const auto& tag_arrs = tags.arrays();

    amrex::ParallelFor(
        ibfab, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const int axp =
                std::abs(ibarrs[nbx](i + 1, j, k) - ibarrs[nbx](i, j, k));
            const int ayp =
                std::abs(ibarrs[nbx](i, j + 1, k) - ibarrs[nbx](i, j, k));
            const int azp =
                std::abs(ibarrs[nbx](i, j, k + 1) - ibarrs[nbx](i, j, k));
            const int axm =
                std::abs(ibarrs[nbx](i - 1, j, k) - ibarrs[nbx](i, j, k));
            const int aym =
                std::abs(ibarrs[nbx](i, j - 1, k) - ibarrs[nbx](i, j, k));
            const int azm =
                std::abs(ibarrs[nbx](i, j, k - 1) - ibarrs[nbx](i, j, k));
            const int ax = amrex::max(axp, axm);
            const int ay = amrex::max(ayp, aym);
            const int az = amrex::max(azp, azm);
            if (amrex::max(ax, ay, az) > 1 ||
                (tag_fringe && ibarrs[nbx](i, j, k) == -1) ||
                (tag_hole && ibarrs[nbx](i, j, k) == 0)) {
                tag_arrs[nbx](i, j, k) = amrex::TagBox::SET;
            }
        });
}

} // namespace amr_wind
