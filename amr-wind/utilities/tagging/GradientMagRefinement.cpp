#include "amr-wind/utilities/tagging/GradientMagRefinement.H"
#include "amr-wind/CFDSim.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

GradientMagRefinement::GradientMagRefinement(const CFDSim& sim)
    : m_sim(sim)
    , m_gradmag_value(
          m_sim.mesh().maxLevel() + 1, std::numeric_limits<amrex::Real>::max())
{}

void GradientMagRefinement::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);
    std::string fname;
    pp.query("field_name", fname);

    const auto& repo = m_sim.repo();
    if (!repo.field_exists(fname)) {
        amrex::Abort("GradientMagRefinement: Cannot find field = " + fname);
    }
    m_field = &(m_sim.repo().get_field(fname));

    amrex::Vector<amrex::Real> gradmag_value;
    pp.queryarr("values", gradmag_value);

    if ((gradmag_value.empty())) {
        amrex::Abort("GradientMagRefinement: Must specify at least one value");
    }

    {
        const int fcount = amrex::min(
            static_cast<int>(gradmag_value.size()),
            static_cast<int>(m_gradmag_value.size()));
        for (int i = 0; i < fcount; ++i) {
            m_gradmag_value[i] = gradmag_value[i];
        }
        m_max_lev_field = fcount - 1;
    }
}

void GradientMagRefinement::operator()(
    int level, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/)
{
    const bool tag_field = level <= m_max_lev_field;
    if (tag_field) {
        m_field->fillpatch(level, time, (*m_field)(level), 1);
    }

    const auto& mfab = (*m_field)(level);
    const auto& idx = m_sim.repo().mesh().Geom(level).InvCellSizeArray();

    const auto& tag_arrs = tags.arrays();
    const auto& farrs = mfab.const_arrays();
    const auto gradmag_val = m_gradmag_value[level];

    amrex::ParallelFor(
        mfab, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            // TODO: ignoring wall stencils for now

            const auto gx =
                0.5_rt * (farrs[nbx](i + 1, j, k) - farrs[nbx](i - 1, j, k)) *
                idx[0];
            const auto gy =
                0.5_rt * (farrs[nbx](i, j + 1, k) - farrs[nbx](i, j - 1, k)) *
                idx[1];
            const auto gz =
                0.5_rt * (farrs[nbx](i, j, k + 1) - farrs[nbx](i, j, k - 1)) *
                idx[2];

            const auto grad_mag = sqrt((gx * gx) + (gy * gy) + (gz * gz));
            if (grad_mag > gradmag_val) {
                tag_arrs[nbx](i, j, k) = amrex::TagBox::SET;
            }
        });
}

} // namespace amr_wind
