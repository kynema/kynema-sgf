#ifndef FILLFLATHER_H
#define FILLFLATHER_H

#include "src/core/FieldFillPatchOps.H"
#include "src/core/FieldBCOps.H"
#include "src/boundary_conditions/field_boundary_fill/Flather.H"

namespace kynema_sgf {

/** Fill patch operator that applies Flather boundary updates. */
class FillFlather : public FieldFillPatchOps<FieldBCDirichlet>
{
public:
    FillFlather(
        Field& field,
        const amrex::AmrCore& mesh,
        const SimTime& time,
        Flather& flather);

    ~FillFlather() override;

    void fillpatch(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        FieldState fstate = FieldState::New) override;

    void fillpatch_from_coarse(
        int /* lev */,
        amrex::Real /* time */,
        amrex::MultiFab& /* mfab */,
        const amrex::IntVect& /* nghost */,
        FieldState /* fstate */) override
    {}

    void fillpatch_sibling_fields(
        int lev,
        amrex::Real time,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& mfabs,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& ffabs,
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& cfabs,
        const amrex::IntVect& nghost,
        const amrex::Vector<amrex::BCRec>& bcrec,
        const amrex::Vector<amrex::BCRec>& /* unused */,
        FieldState fstate = FieldState::New) override;

    void fillphysbc(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        FieldState fstate = FieldState::New) override;

protected:
    Flather& m_flather;
};

} // namespace kynema_sgf

#endif /* FILLFLATHER_H */
