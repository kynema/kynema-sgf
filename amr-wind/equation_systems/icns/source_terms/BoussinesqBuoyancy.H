#ifndef BOUSSINESQBUOYANCY_H
#define BOUSSINESQBUOYANCY_H

#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/transport_models/TransportModel.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde::icns {

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
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const Field& m_temperature;

    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};

    //! Transport model
    const transport::TransportModel& m_transport;
};
} // namespace amr_wind::pde::icns

#endif /* BOUSSINESQBUOYANCY_H */
