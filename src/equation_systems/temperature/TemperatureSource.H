#ifndef TEMPERATURESOURCE_H
#define TEMPERATURESOURCE_H

#include "src/core/Factory.H"
#include "src/core/FieldDescTypes.H"
#include "src/core/FieldUtils.H"
#include "src/core/FieldRepo.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

class CFDSim;

namespace pde {

/** Representation of a temperature source term
 *  \ingroup temp_src
 *
 *  All temperature source terms must inherit from this class
 */
class TemperatureSource : public Factory<TemperatureSource, const CFDSim&>
{
public:
    static std::string base_identifier() { return "TemperatureSource"; }

    ~TemperatureSource() override = default;

    virtual void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const = 0;
};

} // namespace pde
} // namespace amr_wind

#endif /* TEMPERATURESOURCE_H */
