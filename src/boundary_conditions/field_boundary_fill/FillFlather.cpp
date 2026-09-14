#include <utility>

#include "src/boundary_conditions/field_boundary_fill/FillFlather.H"

namespace kynema_sgf {

FillFlather::FillFlather(
    Field& field,
    const amrex::AmrCore& mesh,
    const SimTime& time,
    Flather& flather)
    : FieldFillPatchOps<FieldBCDirichlet>(
          field, mesh, time, FieldInterpolator::CellConsLinear)
    , m_flather(flather)
{}

FillFlather::~FillFlather() = default;

void FillFlather::fillpatch(
    const int lev,
    const amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& /* nghost */,
    const FieldState fstate)
{
    if (m_field.base_name() == "velocity") {
        m_flather.update_flather_variables(lev, fstate);
        m_flather.set_velocity(lev, time, m_field, mfab);
    }
}

void FillFlather::fillphysbc(
    const int lev,
    const amrex::Real time,
    amrex::MultiFab& mfab,
    const amrex::IntVect& /* nghost */,
    const FieldState fstate)
{
    if (m_field.base_name() == "velocity") {
        m_flather.update_flather_variables(lev, fstate);
        m_flather.set_velocity(lev, time, m_field, mfab);
    }
}

void FillFlather::fillpatch_sibling_fields(
    const int lev,
    const amrex::Real time,
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& mfabs,
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& /* ffabs */,
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& /* cfabs */,
    const amrex::IntVect& /* nghost */,
    const amrex::Vector<amrex::BCRec>& /* bcrec */,
    const amrex::Vector<amrex::BCRec>& /* unused */,
    const FieldState fstate)
{
    m_flather.update_flather_variables(lev, fstate, true);
    for (int i = 0; std::cmp_less(i, mfabs.size()); ++i) {
        m_flather.set_velocity(lev, time, m_field, *mfabs[i], 0, i);
    }
}

} // namespace kynema_sgf
