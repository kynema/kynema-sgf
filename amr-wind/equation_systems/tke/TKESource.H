#ifndef TKESOURCE_H
#define TKESOURCE_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/core/FieldDescTypes.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/core/FieldRepo.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

class CFDSim;

namespace pde {

/** Base class for all TKE source terms
 *  \ingroup tke_src
 */
class TKESource : public Factory<TKESource, const CFDSim&>
{
public:
    static std::string base_identifier() { return "TKESource"; }

    ~TKESource() override = default;

    virtual void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const = 0;
};

} // namespace pde
} // namespace amr_wind

#endif /* TKESOURCE_H */
