#ifndef ACTUATORFORCING_H
#define ACTUATORFORCING_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"

namespace amr_wind::pde::icns {

/** Body forces introduced by turbines modeled as actuators in flow field.
 *
 *  \ingroup icns_src actuator
 *
 *  \sa Acuator
 */
class ActuatorForcing : public MomentumSource::Register<ActuatorForcing>
{
public:
    static std::string identifier() { return "ActuatorForcing"; }

    explicit ActuatorForcing(const CFDSim& sim);

    ~ActuatorForcing() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const Field& m_act_src;
};

} // namespace amr_wind::pde::icns

#endif /* ACTUATORFORCING_H */
