#ifndef RAYLEIGHTAYLOR_H
#define RAYLEIGHTAYLOR_H

#include <memory>

#include "amr-wind/core/Physics.H"
#include "amr-wind/physics/RayleighTaylorFieldInit.H"
#include "amr-wind/core/Field.H"

namespace amr_wind {

/** Rayleigh Taylor physics
 *  \ingroup physics
 *
 *  \sa RayleighTaylorFieldInit
 */
class RayleighTaylor : public Physics::Register<RayleighTaylor>
{
public:
    static std::string identifier() { return "RayleighTaylor"; }

    explicit RayleighTaylor(const CFDSim& sim);

    ~RayleighTaylor() override = default;

    //! Initialize the density and velocity fields for RayleighTaylor
    //! simulations
    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void post_init_actions() override {}

    void post_regrid_actions() override {}

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    Field& m_velocity;
    Field& m_density;

    //! RayleighTaylor field initializer instance
    std::unique_ptr<RayleighTaylorFieldInit> m_field_init;
};

} // namespace amr_wind

#endif /* RayleighTaylor_H */
