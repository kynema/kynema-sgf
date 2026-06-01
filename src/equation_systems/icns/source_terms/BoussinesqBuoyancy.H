#ifndef BOUSSINESQBUOYANCY_H
#define BOUSSINESQBUOYANCY_H

#include "src/core/FieldRepo.H"
#include "src/transport_models/TransportModel.H"
#include "src/equation_systems/icns/MomentumSource.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::icns {

/** Boussinesq buoyancy source term
 *  \ingroup icns_src we_abl
 *
 *  \f[
 *    S = \beta g \left( T_\mathrm{ref} - T \right)
 *  \f]
 */
class BoussinesqBuoyancy : public MomentumSource::Register<BoussinesqBuoyancy>
{
public:
    static std::string identifier() { return "BoussinesqBuoyancy"; }

    explicit BoussinesqBuoyancy(const CFDSim& sim);

    ~BoussinesqBuoyancy() override;

    void operator()(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override;

private:
    const Field& m_temperature;
    FieldRepo& m_repo;

    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};

    //! Transport model
    const transport::TransportModel& m_transport;

    //! Scratch fields for level-wide beta and ref_theta
    mutable std::unique_ptr<ScratchField> m_beta_scratch;
    mutable std::unique_ptr<ScratchField> m_ref_theta_scratch;
};
} // namespace kynema_sgf::pde::icns

#endif /* BOUSSINESQBUOYANCY_H */
