#ifndef GRAVITYFORCING_H
#define GRAVITYFORCING_H

#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::icns {

/** Gravity Forcing term
 *  \ingroup icns_src
 */
class GravityForcing : public MomentumSource::Register<GravityForcing>
{
public:
    static std::string identifier() { return "GravityForcing"; }

    explicit GravityForcing(const CFDSim& sim);

    ~GravityForcing() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& vel_forces) const override;

private:
    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};

    Field* m_rho{nullptr};
    Field* m_rho0{nullptr};

    // Input reference to perturbational form
    bool m_use_perturb_pressure{false};
    // Presence of rho0 for perturbational form
    bool m_use_reference_density{false};
    // Constant reference density
    amrex::Real m_rho0_const{1.0_rt};
};

} // namespace amr_wind::pde::icns

#endif /* GravityForcing_H */
