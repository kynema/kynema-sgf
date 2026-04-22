#ifndef INTFIELD_H
#define INTFIELD_H

#include <string>

#include "amr-wind/core/FieldDescTypes.H"

#include "AMReX_iMultiFab.H"

namespace amr_wind {

class FieldRepo;

/** An integer computational field
 *  \ingroup fields
 *
 *  Used to indicate masking and overset state with IBlanking
 */
class IntField
{
public:
    friend class FieldRepo;

    IntField(const IntField&) = delete;
    IntField& operator=(const IntField&) = delete;

    //! Name of the field
    [[nodiscard]] const std::string& name() const { return m_name; }

    //! Unique integer ID for this field
    [[nodiscard]] unsigned id() const { return m_id; }

    //! Number of components for this field
    [[nodiscard]] int num_comp() const { return m_ncomp; }

    //! Number of ghost cells
    [[nodiscard]] const amrex::IntVect& num_grow() const { return m_ngrow; }

    //! Location of the field
    [[nodiscard]] FieldLoc field_location() const { return m_floc; }

    //! Reference to the FieldRepo that holds the fabs
    [[nodiscard]] const FieldRepo& repo() const { return m_repo; }

    //! Access the FAB at a given level
    amrex::iMultiFab& operator()(int lev);
    const amrex::iMultiFab& operator()(int lev) const;

    amrex::Vector<amrex::iMultiFab*> vec_ptrs();

    [[nodiscard]] amrex::Vector<const amrex::iMultiFab*> vec_const_ptrs() const;

    void setVal(int value);

    void setVal(int value, int start_comp, int num_comp = 1, int nghost = 0);

    void setVal(const amrex::Vector<int>& values, int nghost = 0);

protected:
    IntField(
        FieldRepo& repo,
        std::string name,
        unsigned fid,
        int ncomp = 1,
        int ngrow = 1,
        FieldLoc floc = FieldLoc::CELL);

    FieldRepo& m_repo;

    std::string m_name;

    const unsigned m_id;

    int m_ncomp;

    amrex::IntVect m_ngrow;

    FieldLoc m_floc;
};

} // namespace amr_wind

#endif /* INTFIELD_H */
