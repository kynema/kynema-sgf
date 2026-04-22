#ifndef MMSFORCING_H
#define MMSFORCING_H

#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"

namespace amr_wind::pde::icns::mms {

class MMSForcing : public MomentumSource::Register<MMSForcing>
{
public:
    static std::string identifier() { return "MMSForcing"; }

    explicit MMSForcing(const CFDSim& sim);

    ~MMSForcing() override = default;

    void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const Field& m_mms_vel_source;
};
} // namespace amr_wind::pde::icns::mms

#endif /* MMSFORCING_H */
