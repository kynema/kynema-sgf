#ifndef OCEANWAVESFILLINFLOW_H
#define OCEANWAVESFILLINFLOW_H

#include "amr-wind/core/FieldFillPatchOps.H"
#include "amr-wind/core/FieldBCOps.H"
#include "amr-wind/ocean_waves/boundary_ops/OceanWavesBoundary.H"

namespace amr_wind {

class OceanWavesFillInflow : public FieldFillPatchOps<FieldBCDirichlet>
{
public:
    OceanWavesFillInflow(
        Field& field,
        const amrex::AmrCore& mesh,
        const SimTime& time,
        const OceanWavesBoundary& ow_bndry);

    ~OceanWavesFillInflow() override;

    void fillpatch(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        FieldState fstate = FieldState::New) override;

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

    //! Implementation that handles filling patches from a coarse to fine level
    void fillpatch_from_coarse(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        FieldState fstate = FieldState::New) override;

    //! Implementation that handles filling physical boundary conditions
    void fillphysbc(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        FieldState fstate = FieldState::New) override;

protected:
    const OceanWavesBoundary& m_ow_bndry;
};

} // namespace amr_wind

#endif /* ABLFILLINFLOW_H */
