#ifndef SCRATCHFIELD_H
#define SCRATCHFIELD_H

#include <string>
#include <utility>

#include "amr-wind/core/FieldDescTypes.H"
#include "amr-wind/core/ViewField.H"
#include "AMReX_MultiFab.H"
#include "AMReX_Vector.H"
#include "AMReX_PhysBCFunct.H"

namespace amr_wind {

class FieldRepo;

/** A temporary computational field
 *  \ingroup fields
 *
 *  A ScratchField is similar to a Field in usage, but is temporary in nature.
 *  It is used as a scratch buffer to compute intermediate quantities. However,
 *  unlike fields these don't have multiple states, and cannot survive across a
 *  regrid. By default, FieldRepo returns a unique pointer to this instance and
 *  it is not safe to hold this pointer across timesteps.
 *
 *  At present, ScratchField cannot be used for I/O and/or post-processing
 * utilities.
 */
class ScratchField
{
public:
    friend class FieldRepo;

    ScratchField(const ScratchField&) = delete;
    ScratchField& operator=(const ScratchField&) = delete;

    //! Name if available for this scratch field
    [[nodiscard]] const std::string& name() const { return m_name; }

    //! Number of components for this field
    [[nodiscard]] int num_comp() const { return m_ncomp; }

    //! Ghost cells
    [[nodiscard]] const amrex::IntVect& num_grow() const { return m_ngrow; }

    //! Cell, node, face where the field is stored
    [[nodiscard]] FieldLoc field_location() const { return m_floc; }

    //! Return the field data for a given level
    amrex::MultiFab& operator()(int lev) { return m_data[lev]; }
    const amrex::MultiFab& operator()(int lev) const { return m_data[lev]; }

    amrex::Vector<amrex::MultiFab*> vec_ptrs()
    {
        return amrex::GetVecOfPtrs(m_data);
    }

    [[nodiscard]] amrex::Vector<const amrex::MultiFab*> vec_const_ptrs() const
    {
        return amrex::GetVecOfConstPtrs(m_data);
    }

    //! Return a reference to the field repository that created this field
    [[nodiscard]] const FieldRepo& repo() const { return m_repo; }

    /** Return a sub-view of the ScratchField instance
     */
    ViewField<ScratchField> subview(const int scomp = 0, const int ncomp = 1)
    {
        return {*this, scomp, ncomp};
    }

    void fillpatch(amrex::Real time);
    void fillpatch(amrex::Real time, const amrex::IntVect& ng);

    void fillpatch(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost);

    void fillpatch(
        int lev,
        amrex::Real time,
        amrex::MultiFab& mfab,
        const amrex::IntVect& nghost,
        amrex::Vector<amrex::BCRec>& bcrec);

protected:
    ScratchField(
        const FieldRepo& repo,
        std::string name,
        const int ncomp = 1,
        const int ngrow = 1,
        const FieldLoc floc = FieldLoc::CELL)
        : m_repo(repo)
        , m_name(std::move(name))
        , m_ncomp(ncomp)
        , m_ngrow(ngrow)
        , m_floc(floc)
    {}

    const FieldRepo& m_repo;
    std::string m_name;
    int m_ncomp;
    amrex::IntVect m_ngrow;
    FieldLoc m_floc;

    amrex::Vector<amrex::MultiFab> m_data;
};

} // namespace amr_wind

#endif /* SCRATCHFIELD_H */
