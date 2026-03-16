#ifndef SOURCETERM_H
#define SOURCETERM_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/core/FieldDescTypes.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/core/FieldRepo.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {
class CFDSim;

namespace pde {

/** Abstract representation of a source term
 *  \ingroup pdeop
 */
class SourceTerm : public Factory<SourceTerm, const CFDSim&>
{
public:
    static std::string base_identifier() { return "SourceTerm"; }

    ~SourceTerm() override = default;

    virtual void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const = 0;
};

} // namespace pde
} // namespace amr_wind

#endif /* SOURCETERM_H */
