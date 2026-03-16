#include "amr-wind/utilities/tagging/VorticityMagRefinement.H"
#include "amr-wind/CFDSim.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

VorticityMagRefinement::VorticityMagRefinement(const CFDSim& sim)
    : m_sim(sim)
    , m_vort_value(
          m_sim.mesh().maxLevel() + 1, std::numeric_limits<amrex::Real>::max())
{}

void VorticityMagRefinement::initialize(const std::string& key)
{
    std::string fname = "velocity";

    const auto& repo = m_sim.repo();
    if (!repo.field_exists(fname)) {
        amrex::Abort("FieldRefinement: Cannot find field = " + fname);
    }
    m_vel = &(m_sim.repo().get_field(fname));

    amrex::Vector<amrex::Real> vort_value;
    amrex::ParmParse pp(key);

    pp.queryarr("values", vort_value);

    if (vort_value.empty()) {
        amrex::Abort(
            "VorticityMagRefinement: Must specify at least one of value");
    }

    {
        const int fcount = amrex::min(
            static_cast<int>(vort_value.size()),
            static_cast<int>(m_vort_value.size()));
        for (int i = 0; i < fcount; ++i) {
            m_vort_value[i] = vort_value[i];
        }
        m_max_lev_field = fcount - 1;
    }
}

void VorticityMagRefinement::operator()(
    int level, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/)
{
    const bool tag_field = level <= m_max_lev_field;

    if (!tag_field) {
        return;
    }

    m_vel->fillpatch(level, time, (*m_vel)(level), 1);

    const auto& mfab = (*m_vel)(level);
    const auto& idx = m_sim.repo().mesh().Geom(level).InvCellSizeArray();

    const auto& tag_arrs = tags.arrays();
    const auto& vel_arrs = mfab.const_arrays();
    const auto vort_val = m_vort_value[level];

    amrex::ParallelFor(
        mfab, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            // TODO: ignoring wall stencils for now
            const auto vx = 0.5_rt *
                            (vel_arrs[nbx](i + 1, j, k, 1) -
                             vel_arrs[nbx](i - 1, j, k, 1)) *
                            idx[0];
            const auto wx = 0.5_rt *
                            (vel_arrs[nbx](i + 1, j, k, 2) -
                             vel_arrs[nbx](i - 1, j, k, 2)) *
                            idx[0];

            const auto uy = 0.5_rt *
                            (vel_arrs[nbx](i, j + 1, k, 0) -
                             vel_arrs[nbx](i, j - 1, k, 0)) *
                            idx[1];
            const auto wy = 0.5_rt *
                            (vel_arrs[nbx](i, j + 1, k, 2) -
                             vel_arrs[nbx](i, j - 1, k, 2)) *
                            idx[1];

            const auto uz = 0.5_rt *
                            (vel_arrs[nbx](i, j, k + 1, 0) -
                             vel_arrs[nbx](i, j, k - 1, 0)) *
                            idx[2];
            const auto vz = 0.5_rt *
                            (vel_arrs[nbx](i, j, k + 1, 1) -
                             vel_arrs[nbx](i, j, k - 1, 1)) *
                            idx[2];

            const auto vort = sqrt(
                std::pow(uy - vx, 2) + std::pow(vz - wy, 2) +
                std::pow(wx - uz, 2));

            if (vort > vort_val) {
                tag_arrs[nbx](i, j, k) = amrex::TagBox::SET;
            }
        });
}

} // namespace amr_wind
