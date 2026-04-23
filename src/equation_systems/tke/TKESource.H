#ifndef TKESOURCE_H
#define TKESOURCE_H

#include "src/core/Factory.H"
#include "src/core/FieldDescTypes.H"
#include "src/core/FieldUtils.H"
#include "src/core/FieldRepo.H"
#include "AMReX_MultiFab.H"

namespace kynema_sgf {

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
} // namespace kynema_sgf

#endif /* TKESOURCE_H */
