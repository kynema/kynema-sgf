#ifndef SYNTHTURBFORCING_H
#define SYNTHTURBFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"

namespace amr_wind::pde::icns {

/** Synthetic turbulence momentum source term
 *  \ingroup icns_src synth_turb
 */
class SynthTurbForcing : public MomentumSource::Register<SynthTurbForcing>
{
public:
    static std::string identifier() { return "SynthTurbForcing"; }

    explicit SynthTurbForcing(const CFDSim& /*sim*/);

    ~SynthTurbForcing() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    Field& m_turb_force;
};

} // namespace amr_wind::pde::icns

#endif /* SYNTHTURBFORCING_H */
