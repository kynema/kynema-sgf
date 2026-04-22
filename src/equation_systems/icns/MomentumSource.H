#ifndef MOMENTUMSOURCE_H
#define MOMENTUMSOURCE_H

#include "src/core/Factory.H"
#include "src/core/FieldDescTypes.H"
#include "src/core/FieldUtils.H"
#include "src/core/FieldRepo.H"
#include "AMReX_MultiFab.H"

namespace kynema_sgf {

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
} // namespace kynema_sgf

#endif /* MOMENTUMSOURCE_H */
