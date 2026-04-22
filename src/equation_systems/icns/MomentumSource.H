#ifndef MOMENTUMSOURCE_H
#define MOMENTUMSOURCE_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/core/FieldDescTypes.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/core/FieldRepo.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

class CFDSim;

namespace pde {

/** Representation of a momentum source term
 *  \ingroup icns_src
 *
 *  All momentum source terms must inherit from this class
 */
class MomentumSource : public Factory<MomentumSource, const CFDSim&>
{
public:
    static std::string base_identifier() { return "MomentumSource"; }

    ~MomentumSource() override = default;

    virtual void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const = 0;
};

} // namespace pde
} // namespace amr_wind

#endif /* MOMENTUMSOURCE_H */
