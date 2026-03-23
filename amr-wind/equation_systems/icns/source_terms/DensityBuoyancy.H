#ifndef DENSITYBUOYANCY_H
#define DENSITYBUOYANCY_H

#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::icns {

/** Buoyancy term from density variations
 *  \ingroup icns_src
 *
 *  \sa RayleighTaylor
 */
class DensityBuoyancy : public MomentumSource::Register<DensityBuoyancy>
{
public:
    static std::string identifier() { return "DensityBuoyancy"; }

    explicit DensityBuoyancy(const CFDSim& sim);

    ~DensityBuoyancy() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& vel_forces) const override;

private:
    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};

    amrex::Real m_rho_0{1.0_rt};

    Field& m_density;
};

} // namespace amr_wind::pde::icns

#endif /* DensityBuoyancy_H */
