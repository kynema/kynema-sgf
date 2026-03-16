#include "amr-wind/utilities/tagging/FieldRefinement.H"
#include "amr-wind/CFDSim.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

FieldRefinement::FieldRefinement(const CFDSim& sim)
    : m_sim(sim)
    , m_field_error(
          m_sim.mesh().maxLevel() + 1, std::numeric_limits<amrex::Real>::max())
    , m_grad_error(
          m_sim.mesh().maxLevel() + 1, std::numeric_limits<amrex::Real>::max())
    , m_tagging_box(m_sim.repo().mesh().Geom(0).ProbDomain())
{}

void FieldRefinement::initialize(const std::string& key)
{
    amrex::ParmParse pp(key);
    std::string fname;
    pp.query("field_name", fname);

    const auto& repo = m_sim.repo();
    if (repo.field_exists(fname)) {
        m_field = &(m_sim.repo().get_field(fname));
        AMREX_ALWAYS_ASSERT(m_field->num_grow() > amrex::IntVect{0});
    } else if (repo.int_field_exists(fname)) {
        m_int_field = &(m_sim.repo().get_int_field(fname));
        AMREX_ALWAYS_ASSERT(m_int_field->num_grow() > amrex::IntVect{0});
    } else {
        amrex::Abort("FieldRefinement: Cannot find field = " + fname);
    }

    amrex::Vector<amrex::Real> field_err;
    amrex::Vector<amrex::Real> grad_err;
    pp.queryarr("field_error", field_err);
    pp.queryarr("grad_error", grad_err);

    amrex::Vector<amrex::Real> box_lo(AMREX_SPACEDIM, 0);
    amrex::Vector<amrex::Real> box_hi(AMREX_SPACEDIM, 0);
    if (pp.queryarr("box_lo", box_lo, 0, static_cast<int>(box_lo.size())) ==
        1) {
        m_tagging_box.setLo(box_lo);
    }
    if (pp.queryarr("box_hi", box_hi, 0, static_cast<int>(box_hi.size())) ==
        1) {
        m_tagging_box.setHi(box_hi);
    }

    if ((field_err.empty()) && (grad_err.empty())) {
        amrex::Abort(
            "FieldRefinement: Must specify at least one of field_error or "
            "grad_error");
    }

    {
        const int fcount = amrex::min(
            static_cast<int>(field_err.size()),
            static_cast<int>(m_field_error.size()));
        for (int i = 0; i < fcount; ++i) {
            m_field_error[i] = field_err[i];
        }
        m_max_lev_field = fcount - 1;
    }
    {
        const int fcount = amrex::min(
            static_cast<int>(grad_err.size()),
            static_cast<int>(m_grad_error.size()));
        for (int i = 0; i < fcount; ++i) {
            m_grad_error[i] = grad_err[i];
        }
        m_max_lev_grad = fcount - 1;
    }
}

void FieldRefinement::operator()(
    const int level,
    amrex::TagBoxArray& tags,
    const amrex::Real time,
    const int /*ngrow*/)
{
    const bool tag_grad = level <= m_max_lev_grad;
    if (tag_grad) {
        if (m_field != nullptr) {
            m_field->fillpatch(level, time, (*m_field)(level), 1);
        } else if (m_int_field != nullptr) {
            (*m_int_field)(level).FillBoundary(
                m_sim.repo().mesh().Geom(level).periodicity());
        }
    }

    if (m_field != nullptr) {
        const auto& mfab = (*m_field)(level);
        tag(level, tags, mfab);
    } else if (m_int_field != nullptr) {
        const auto& mfab = (*m_int_field)(level);
        tag(level, tags, mfab);
    }
}

} // namespace amr_wind
