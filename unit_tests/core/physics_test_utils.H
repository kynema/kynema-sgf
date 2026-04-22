#ifndef PHYSICS_TEST_UTILS_H
#define PHYSICS_TEST_UTILS_H

#include <string>

#include "src/core/Physics.H"

namespace amr_wind_tests {

class PhysicsEx : public amr_wind::Physics::Register<PhysicsEx>
{
public:
    static std::string identifier() { return "PhysicsEx"; }

    explicit PhysicsEx(amr_wind::CFDSim& /*unused*/);

    ~PhysicsEx() override = default;

    void post_init_actions() override;
    void post_regrid_actions() override;
    void
    initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/) override;
    void pre_advance_work() override;
    void post_advance_work() override;
};

} // namespace amr_wind_tests

#endif /* PHYSICS_TEST_UTILS_H */
