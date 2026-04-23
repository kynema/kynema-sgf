#ifndef BURGGRAF_FLOW_FORCING_H
#define BURGGRAF_FLOW_FORCING_H

#include "src/equation_systems/icns/MomentumSource.H"
#include "src/core/SimTime.H"

namespace kynema_sgf::pde::icns {

/** Body forcing introduced to represent the Burrgraf cavity flow.
 *
 *  \ingroup icns_src burggraf_flow
 */
class BurggrafFlowForcing : public MomentumSource::Register<BurggrafFlowForcing>
{
public:
    static std::string identifier() { return "BurggrafFlowForcing"; }

    explicit BurggrafFlowForcing(const CFDSim& sim);

    ~BurggrafFlowForcing() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const Field& m_bf_src;
};

} // namespace kynema_sgf::pde::icns

#endif /* ACTUATORFORCING_H */
