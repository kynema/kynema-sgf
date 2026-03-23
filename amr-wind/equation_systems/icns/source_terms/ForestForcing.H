#ifndef FORESTFORCING_H
#define FORESTFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::pde::icns {

/** Adds a forcing term to represent forest terrain
 *
 *  \ingroup icns_src
 *
 *
 */
class ForestForcing : public MomentumSource::Register<ForestForcing>
{
public:
    static std::string identifier() { return "ForestForcing"; }

    explicit ForestForcing(const CFDSim& sim);

    ~ForestForcing() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const CFDSim& m_sim;
    const Field& m_velocity;
};

} // namespace amr_wind::pde::icns

#endif
