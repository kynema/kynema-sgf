#include "mms_test_utils.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"

#include "amr-wind/physics/mms/MMS.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {
TEST_F(MMSMeshTest, mms_initialization)
{
#if defined(AMREX_USE_HIP)
    GTEST_SKIP();
#else
    populate_parameters();
    utils::populate_mms_params();

    initialize_mesh();
    auto& frepo = mesh().field_repo();
    auto& velocityf = frepo.declare_field("velocity", 3, 0);
    auto& densityf = frepo.declare_field("density");

    amr_wind::mms::MMS mms(sim());
    const int nlevels = mesh().num_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        mms.initialize_fields(lev, mesh().Geom(lev));
    }

    auto velocity = velocityf.vec_ptrs();
    auto density = densityf.vec_ptrs();
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    // Test density
    {
        amrex::Real min_rho, max_rho;
        utils::field_minmax(nlevels, density, min_rho, max_rho);
        EXPECT_NEAR(min_rho, 1.0_rt, tol);
        EXPECT_NEAR(max_rho, 1.0_rt, tol);
    }

    // Test velocity
    {
        amrex::Vector<amrex::Real> min_vel(3), max_vel(3);
        utils::field_minmax(nlevels, velocity, min_vel, max_vel);
        EXPECT_NEAR(min_vel[0], -0.085642667425670935_rt, tol);
        EXPECT_NEAR(min_vel[1], -0.15852969775284545_rt, tol);
        EXPECT_NEAR(min_vel[2], -0.28208275785136977_rt, tol);
        EXPECT_NEAR(max_vel[0], 0.085642667425670921_rt, tol);
        EXPECT_NEAR(max_vel[1], 0.15852969775284548_rt, tol);
        EXPECT_NEAR(max_vel[2], 0.28208275785136988_rt, tol);
    }
#endif
}

} // namespace amr_wind_tests
