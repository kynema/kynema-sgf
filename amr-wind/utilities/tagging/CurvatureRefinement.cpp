#include "amr-wind/utilities/tagging/CurvatureRefinement.H"
#include "amr-wind/CFDSim.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

CurvatureRefinement::CurvatureRefinement(const CFDSim& sim)
    : m_sim(sim)
    , m_curv_value(
          m_sim.mesh().maxLevel() + 1, std::numeric_limits<amrex::Real>::max())
{}

void CurvatureRefinement::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);
    std::string fname;
    pp.query("field_name", fname);

    const auto& repo = m_sim.repo();
    if (!repo.field_exists(fname)) {
        amrex::Abort("CurvatureRefinement: Cannot find field = " + fname);
    }
    m_field = &(m_sim.repo().get_field(fname));

    amrex::Vector<amrex::Real> curv_value;
    pp.queryarr("values", curv_value);

    if ((curv_value.empty())) {
        amrex::Abort("CurvatureRefinement: Must specify at least one value");
    }

    {
        const int fcount = amrex::min(
            static_cast<int>(curv_value.size()),
            static_cast<int>(m_curv_value.size()));
        for (int i = 0; i < fcount; ++i) {
            m_curv_value[i] = curv_value[i];
        }
        m_max_lev_field = fcount - 1;
    }
}

void CurvatureRefinement::operator()(
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
    const auto curv_val = m_curv_value[level];

    amrex::ParallelFor(
        mfab, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            // TODO: ignoring wall stencils for now

            const auto phixx =
                (farrs[nbx](i + 1, j, k) - 2.0_rt * farrs[nbx](i, j, k) +
                 farrs[nbx](i - 1, j, k)) *
                idx[0] * idx[0];
            const auto phiyy =
                (farrs[nbx](i, j + 1, k) - 2.0_rt * farrs[nbx](i, j, k) +
                 farrs[nbx](i, j - 1, k)) *
                idx[0] * idx[0];
            const auto phizz =
                (farrs[nbx](i, j, k + 1) - 2.0_rt * farrs[nbx](i, j, k) +
                 farrs[nbx](i, j, k - 1)) *
                idx[0] * idx[0];

            const auto phiz =
                0.5_rt * (farrs[nbx](i, j, k + 1) - farrs[nbx](i, j, k - 1)) *
                idx[2];
            const auto phiz_ip1 =
                0.5_rt *
                (farrs[nbx](i + 1, j, k + 1) - farrs[nbx](i + 1, j, k - 1)) *
                idx[2];
            const auto phiz_im1 =
                0.5_rt *
                (farrs[nbx](i - 1, j, k + 1) - farrs[nbx](i - 1, j, k - 1)) *
                idx[2];
            const auto phiz_jp1 =
                0.5_rt *
                (farrs[nbx](i, j + 1, k + 1) - farrs[nbx](i, j + 1, k - 1)) *
                idx[2];
            const auto phiz_jm1 =
                0.5_rt *
                (farrs[nbx](i, j - 1, k + 1) - farrs[nbx](i, j - 1, k - 1)) *
                idx[2];

            const auto phiy =
                0.5_rt * (farrs[nbx](i, j + 1, k) - farrs[nbx](i, j - 1, k)) *
                idx[1];
            const auto phiy_ip1 =
                0.5_rt *
                (farrs[nbx](i + 1, j + 1, k) - farrs[nbx](i + 1, j - 1, k)) *
                idx[1];
            const auto phiy_im1 =
                0.5_rt *
                (farrs[nbx](i - 1, j + 1, k) - farrs[nbx](i - 1, j - 1, k)) *
                idx[1];
            const auto phiyz = 0.5_rt * (phiz_jp1 - phiz_jm1) * idx[1];

            const auto phix =
                0.5_rt * (farrs[nbx](i + 1, j, k) - farrs[nbx](i - 1, j, k)) *
                idx[0];
            const auto phixy = 0.5_rt * (phiy_ip1 - phiy_im1) * idx[0];
            const auto phixz = 0.5_rt * (phiz_ip1 - phiz_im1) * idx[0];

            const auto curv_mag =
                std::abs(
                    (phix * phix * phiyy) - (2. * phix * phiy * phixy) +
                    (phiy * phiy * phixx) + (phix * phix * phizz) -
                    (2. * phix * phiz * phixz) + (phiz * phiz * phixx) +
                    (phiy * phiy * phizz) - (2. * phiy * phiz * phiyz) +
                    (phiz * phiz * phiyy)) /
                std::pow((phix * phix) + (phiy * phiy) + (phiz * phiz), 1.5_rt);
            const auto curv_min = amrex::min<amrex::Real>(
                curv_val, std::cbrt(idx[0] * idx[1] * idx[2]));
            if (curv_mag > curv_min) {
                tag_arrs[nbx](i, j, k) = amrex::TagBox::SET;
            }
        });
}

} // namespace amr_wind
