#ifndef ABLFILLINFLOW_H
#define ABLFILLINFLOW_H

#include "amr-wind/core/FieldFillPatchOps.H"
#include "amr-wind/core/FieldBCOps.H"
#include "amr-wind/wind_energy/ABLBoundaryPlane.H"

namespace amr_wind {

/** Operations to perform inflow of data from boundary planes
 *  \ingroup we_abl
 *
 *  This class contains the operations necessary to interface with
 *  ABLBoundaryPlane class and populate data at the boundaries
 *
 *  \sa ABLBoundaryPlane
 */
class ABLFillInflow : public FieldFillPatchOps<FieldBCDirichlet>
{
public:
    ABLFillInflow(
        Field& field,
        const amrex::AmrCore& mesh,
        const SimTime& time,
        const ABLBoundaryPlane& bndry_plane);

    ~ABLFillInflow() override;

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
    const ABLBoundaryPlane& m_bndry_plane;
};

} // namespace amr_wind

#endif /* ABLFILLINFLOW_H */
