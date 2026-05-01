#include "physics_test_utils.H"

namespace kynema_sgf_tests {

PhysicsEx::PhysicsEx(kynema_sgf::CFDSim& /*unused*/) {}

void PhysicsEx::post_init_actions() {}
void PhysicsEx::post_regrid_actions() {}
void PhysicsEx::initialize_fields(
    int /*level*/, const amrex::Geometry& /*geom*/)
{}
void PhysicsEx::pre_advance_work() {}
void PhysicsEx::post_advance_work() {}

} // namespace kynema_sgf_tests
