#ifndef PHYSICS_TEST_UTILS_H
#define PHYSICS_TEST_UTILS_H

#include <string>

#include "src/core/Physics.H"

namespace kynema_sgf_tests {

class PhysicsEx : public kynema_sgf::Physics::Register<PhysicsEx>
{
public:
    static std::string identifier() { return "PhysicsEx"; }

    explicit PhysicsEx(kynema_sgf::CFDSim& /*unused*/);

    ~PhysicsEx() override = default;

    void post_init_actions() override;
    void post_regrid_actions() override;
    void
    initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/) override;
    void pre_advance_work() override;
    void post_advance_work() override;
};

} // namespace kynema_sgf_tests

#endif /* PHYSICS_TEST_UTILS_H */
