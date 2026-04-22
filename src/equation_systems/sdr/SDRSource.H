#ifndef SDRSOURCE_H
#define SDRSOURCE_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/core/FieldDescTypes.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/core/FieldRepo.H"
#include "AMReX_MultiFab.H"

namespace amr_wind {

class CFDSim;

namespace pde {

/** Base class for SDR source terms.
 *  \ingroup sdr_src
 *
 *  All concrete implementations of SDR source terms must subclass this class.
 */
class SDRSource : public Factory<SDRSource, const CFDSim&>
{
public:
    static std::string base_identifier() { return "SDRSource"; }

    ~SDRSource() override = default;

    virtual void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const = 0;
};

} // namespace pde
} // namespace amr_wind

#endif /* SDRSOURCE_H */
