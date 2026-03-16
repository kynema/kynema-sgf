#ifndef ABLFILLMPL_H
#define ABLFILLMPL_H

#include "amr-wind/core/FieldFillPatchOps.H"
#include "amr-wind/core/FieldBCOps.H"
#include "amr-wind/wind_energy/ABLModulatedPowerLaw.H"

namespace amr_wind {

/** Operations to perform inflow of data from boundary planes
 *  \ingroup we_abl
 *
 *  This class contains the operations necessary to interface with
 *  ABLBoundaryPlane class and populate data at the boundaries
 *
 *  \sa ABLBoundaryPlane
 */
class ABLFillMPL : public FieldFillPatchOps<FieldBCDirichlet>
{
public:
    ABLFillMPL(
        Field& field,
        const amrex::AmrCore& mesh,
        const SimTime& time,
        const ABLModulatedPowerLaw& abl_mpl);

    ~ABLFillMPL() override;

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
    const ABLModulatedPowerLaw& m_abl_mpl;
};

} // namespace amr_wind

#endif /* ABLFILLMPL_H */
