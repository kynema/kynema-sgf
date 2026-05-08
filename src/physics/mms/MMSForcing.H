#ifndef MMSFORCING_H
#define MMSFORCING_H

#include "src/core/FieldRepo.H"
#include "src/equation_systems/icns/MomentumSource.H"
#include "AMReX_MultiFab.H"

namespace kynema_sgf::pde::icns::mms {

class MMSForcing : public MomentumSource::Register<MMSForcing>
{
public:
    static std::string identifier() { return "MMSForcing"; }

    explicit MMSForcing(const CFDSim& sim);

    ~MMSForcing() override = default;

    void operator()(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override;

private:
    const Field& m_mms_vel_source;
};
} // namespace kynema_sgf::pde::icns::mms

#endif /* MMSFORCING_H */
