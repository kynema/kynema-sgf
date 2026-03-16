#include "abl_test_utils.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/wind_energy/ABLFieldInit.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {
namespace {} // namespace

TEST_F(ABLMeshTest, abl_initialization)
{
    populate_parameters();

    // For the first test disable velocity perturbations
    {
        amrex::ParmParse pp("ABL");
        bool perturb_velocity = false;
        pp.add("perturb_velocity", static_cast<int>(perturb_velocity));
    }

    initialize_mesh();
    auto& frepo = mesh().field_repo();
    auto& velocityf = frepo.declare_field("velocity", 3, 0);
    auto& densityf = frepo.declare_field("density");
    auto& temperaturef = frepo.declare_field("temperature");

    auto velocity = velocityf.vec_ptrs();
    auto density = densityf.vec_ptrs();
    auto temperature = temperaturef.vec_ptrs();

    amr_wind::ABLFieldInit ablinit;
    run_algorithm(
        mesh().num_levels(), density,
        [&](const int lev, const amrex::MFIter& mfi) {
            auto vel = velocity[lev]->array(mfi);
            auto rho = density[lev]->array(mfi);
            auto theta = temperature[lev]->array(mfi);

            const auto& bx = mfi.validbox();
            ablinit(bx, mesh().Geom(lev), vel, rho, theta);
        });

    const int nlevels = mesh().num_levels();
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    // Test temperature
    {
        const amrex::Real dz = mesh().Geom(0).CellSize(2);
        const amrex::Real min_temp_gold = 300.0_rt;
        const amrex::Real max_temp_gold =
            308.75_rt - ((308.75_rt - 308.0_rt) / 250.0_rt * (0.5_rt * dz));
        amrex::Real min_temp, max_temp;
        utils::field_minmax(nlevels, temperature, min_temp, max_temp);
        EXPECT_NEAR(min_temp, min_temp_gold, tol);
        EXPECT_NEAR(max_temp, max_temp_gold, tol);
    }

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
        EXPECT_NEAR(min_vel[0], 20.0_rt, tol);
        EXPECT_NEAR(min_vel[1], 10.0_rt, tol);
        EXPECT_NEAR(min_vel[2], 0.0_rt, tol);
        EXPECT_NEAR(max_vel[0], 20.0_rt, tol);
        EXPECT_NEAR(max_vel[1], 10.0_rt, tol);
        EXPECT_NEAR(max_vel[2], 0.0_rt, tol);
    }
}

} // namespace amr_wind_tests
