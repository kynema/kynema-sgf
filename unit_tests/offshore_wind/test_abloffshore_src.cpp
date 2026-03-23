#include <numbers>
#include "abloffshore_test_utils.H"
#include "amr-wind/utilities/trig_ops.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "AMReX_REAL.H"
#include "amr-wind/equation_systems/icns/icns.H"
#include "amr-wind/equation_systems/icns/icns_ops.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/equation_systems/icns/source_terms/ABLForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/GeostrophicForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/BoussinesqBuoyancy.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"

using namespace amrex::literals;

namespace amr_wind_tests {

namespace {
amrex::Real get_val_at_height(
    amr_wind::Field& field,
    const int lev,
    const int comp,
    const amrex::Real ploz,
    const amrex::Real dz,
    const amrex::Real height)
{
    amrex::Real error_total = 0;

    error_total += amrex::ReduceSum(
        field(lev), 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx,
            amrex::Array4<amrex::Real const> const& f_arr) -> amrex::Real {
            amrex::Real error = 0;

            amrex::Loop(bx, [=, &error](int i, int j, int k) {
                const amrex::Real z = ploz + ((0.5_rt + k) * dz);
                // Check if current cell is closest to desired height
                if (z - height < 0.5_rt * dz && z - height >= -0.5_rt * dz) {
                    // Add field value to output
                    error += f_arr(i, j, k, comp);
                }
            });

            return error;
        });
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    return error_total;
}
void init_abl_temperature_field(
    const amrex::Geometry& geom,
    const amrex::Box& bx,
    const amrex::Array4<amrex::Real>& trac,
    const amrex::Real bottom)
{
    const auto& dx = geom.CellSizeArray();
    const auto& plo = geom.ProbLoArray();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real z = plo[2] + ((k + 0.5_rt) * dx[2]);

        // potential temperature profile
        if (z < bottom) {
            trac(i, j, k, 0) = 300.0_rt;
        } else if (z < 250.0_rt) {
            trac(i, j, k, 0) =
                300.0_rt + ((z - bottom) / (250.0_rt - bottom) * 8.0_rt);
        } else {
            trac(i, j, k, 0) = 308.0_rt;
        }
    });
}
void init_vof_field(
    const amrex::Geometry& geom,
    const amrex::Box& bx,
    const amrex::Array4<amrex::Real>& vof,
    const amrex::Real wlev)
{
    const auto& dx = geom.CellSizeArray();
    const auto& plo = geom.ProbLoArray();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real z_btm = plo[2] + (k * dx[2]);
        const amrex::Real z_top = plo[2] + ((k + 1) * dx[2]);

        // vof profile for flat interface
        if (z_btm >= wlev) {
            vof(i, j, k) = 0.0_rt;
        } else if (z_top <= wlev) {
            vof(i, j, k) = 1.0_rt;
        } else {
            vof(i, j, k) = (wlev - z_btm) / (z_top - z_btm);
        }
    });
}
} // namespace

using ICNSFields =
    amr_wind::pde::FieldRegOp<amr_wind::pde::ICNS, amr_wind::fvm::Godunov>;

TEST_F(ABLOffshoreMeshTest, abl_forcing)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    populate_parameters();
    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    pde_mgr.register_transport_pde("Temperature");
    sim().init_physics();
    auto& mphase = sim().physics_manager().get<amr_wind::MultiPhase>();
    // Make sure to read water level
    mphase.post_init_actions();

    auto& src_term = pde_mgr.icns().fields().src_term;

    amr_wind::pde::icns::ABLForcing abl_forcing(sim());

    src_term.setVal(0.0_rt);
    // Mimic source term at later timesteps
    {
        auto& time = sim().time();
        time.new_timestep();
        time.set_current_cfl(2.0_rt, 0.0_rt, 0.0_rt);
        EXPECT_NEAR(time.delta_t(), 0.1_rt, tol);

        src_term.setVal(0.0_rt);
        abl_forcing.set_mean_velocities(10.0_rt, 5.0_rt);
        auto& volume_fraction = sim().repo().get_field("vof");
        auto waterlev = mphase.water_level();
        auto& geom = sim().mesh().Geom();
        // Initialize vof field before calculating forcing function
        run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
            const auto& bx = mfi.tilebox();
            const auto& gbx = mfi.growntilebox(3);
            const auto& vof_arr = volume_fraction(lev).array(mfi);
            init_vof_field(geom[lev], gbx, vof_arr, waterlev);
            const auto& src_arr = src_term(lev).array(mfi);
            abl_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
        });

        // Targets are U = (20.0_rt, 10.0_rt, 0.0_rt) set in initial conditions
        // Means (set above) V = (10.0_rt, 5.0_rt, 0.0_rt)
        // delta_t (set above) dt = 0.1_rt
        const amrex::Array<amrex::Real, AMREX_SPACEDIM> golds{
            {100.0_rt, 50.0_rt, 0.0_rt}};
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            const auto min_val = utils::field_min(src_term, i);
            const auto max_val = utils::field_max(src_term, i);
            // Ensure that the source term value is present
            EXPECT_NEAR(max_val, golds[i], tol);
            // Ensure that the source term is turned off somewhere
            EXPECT_NEAR(min_val, 0.0_rt, tol);
        }

        // Check that values are turned off in off region
        // Check that values are turned off in proximity to interface
        // Check that values follow cosine curve in ramp region
        const amrex::Real dz = geom[0].CellSize(2);
        const amrex::Real ploz = geom[0].ProbLo(2);
        // Forcing limit parameters from abloffshore_test_utils
        const amrex::Real ht0 = 10.0_rt;
        const amrex::Real ht1 = 30.0_rt;
        const amrex::Array<amrex::Real, 5> test_heights{
            waterlev - dz, waterlev + (0.5_rt * ht0), waterlev + (2.0_rt * dz),
            (std::floor((waterlev + ht0 + ht1 - dz) / dz) * dz) + (0.5_rt * dz),
            waterlev + ht0 + ht1 + dz};
        // Expected values of coeff for each location
        const amrex::Array<amrex::Real, 5> coeff_golds{
            0.0_rt, 0.0_rt, 0.0_rt,
            (-0.5_rt * std::cos(
                           std::numbers::pi_v<amrex::Real> *
                           (test_heights[3] - waterlev - ht0) / ht1)) +
                0.5_rt,
            1.0_rt};

        // Get values from src term
        amrex::Array<amrex::Real, 5> src_x_vals;
        amrex::Array<amrex::Real, 5> src_y_vals;
        amrex::Array<amrex::Real, 5> src_z_vals;
        // Divide by nx*ny cells because of sum
        const int nx = 8;
        const int ny = 8;
        for (int n = 0; n < 5; ++n) {
            src_x_vals[n] =
                get_val_at_height(src_term, 0, 0, ploz, dz, test_heights[n]) /
                nx / ny;
            src_y_vals[n] =
                get_val_at_height(src_term, 0, 1, ploz, dz, test_heights[n]) /
                nx / ny;
            src_z_vals[n] =
                get_val_at_height(src_term, 0, 2, ploz, dz, test_heights[n]) /
                nx / ny;
        }

        // Check each src value against expectations
        for (int n = 0; n < 5; ++n) {
            EXPECT_NEAR(src_x_vals[n], coeff_golds[n] * golds[0], tol);
            EXPECT_NEAR(src_y_vals[n], coeff_golds[n] * golds[1], tol);
            EXPECT_NEAR(src_z_vals[n], coeff_golds[n] * golds[2], tol);
        }
    }
}

TEST_F(ABLOffshoreMeshTest, geostrophic_forcing)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    populate_parameters();

    amrex::ParmParse pp("CoriolisForcing");
    pp.add("latitude", 54.0_rt);

    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();
    auto& mphase = sim().physics_manager().get<amr_wind::MultiPhase>();
    // Make sure to read water level
    mphase.post_init_actions();

    auto& src_term = pde_mgr.icns().fields().src_term;
    auto& density = sim().repo().get_field("density");
    density.setVal(1.0_rt);

    amr_wind::pde::icns::GeostrophicForcing geostrophic_forcing(sim());
    src_term.setVal(0.0_rt);
    auto& volume_fraction = sim().repo().get_field("vof");
    auto waterlev = mphase.water_level();
    auto& geom = sim().mesh().Geom();
    run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& gbx = mfi.growntilebox(3);
        const auto& vof_arr = volume_fraction(lev).array(mfi);
        init_vof_field(geom[lev], gbx, vof_arr, waterlev);
        const auto& src_arr = src_term(lev).array(mfi);
        geostrophic_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });

    constexpr amrex::Real corfac =
        2.0_rt * amr_wind::utils::two_pi() / 86164.091_rt * 0.80901699437_rt;
    const amrex::Array<amrex::Real, AMREX_SPACEDIM> golds{
        {-corfac * 6.0_rt, corfac * 10.0_rt, 0.0_rt}};
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const auto min_val = utils::field_min(src_term, i);
        const auto max_val = utils::field_max(src_term, i);
        EXPECT_NEAR(
            amrex::max<amrex::Real>(std::abs(max_val), std::abs(min_val)),
            std::abs(golds[i]), tol);
        EXPECT_NEAR(
            amrex::min<amrex::Real>(std::abs(max_val), std::abs(min_val)),
            0.0_rt, tol);
    }

    // Check that values are turned off in off region
    // Check that values are turned off in proximity to interface
    // Check that values follow cosine curve in ramp region
    const amrex::Real dz = geom[0].CellSize(2);
    const amrex::Real ploz = geom[0].ProbLo(2);
    // Forcing limit parameters from abloffshore_test_utils
    const amrex::Real ht0 = 10.0_rt;
    const amrex::Real ht1 = 30.0_rt;
    const amrex::Array<amrex::Real, 5> test_heights{
        waterlev - dz, waterlev + (0.5_rt * ht0), waterlev + (2.0_rt * dz),
        (std::floor((waterlev + ht0 + ht1 - dz) / dz) * dz) + (0.5_rt * dz),
        waterlev + ht0 + ht1 + dz};
    // Expected values of coeff for each location
    const amrex::Array<amrex::Real, 5> coeff_golds{
        0.0_rt, 0.0_rt, 0.0_rt,
        (-0.5_rt * std::cos(
                       std::numbers::pi_v<amrex::Real> *
                       (test_heights[3] - waterlev - ht0) / ht1)) +
            0.5_rt,
        1.0_rt};

    // Get values from src term
    amrex::Array<amrex::Real, 5> src_x_vals;
    amrex::Array<amrex::Real, 5> src_y_vals;
    amrex::Array<amrex::Real, 5> src_z_vals;
    // Divide by nx*ny cells because of sum
    const int nx = 8;
    const int ny = 8;
    for (int n = 0; n < 5; ++n) {
        src_x_vals[n] =
            get_val_at_height(src_term, 0, 0, ploz, dz, test_heights[n]) / nx /
            ny;
        src_y_vals[n] =
            get_val_at_height(src_term, 0, 1, ploz, dz, test_heights[n]) / nx /
            ny;
        src_z_vals[n] =
            get_val_at_height(src_term, 0, 2, ploz, dz, test_heights[n]) / nx /
            ny;
    }

    // Check each src value against expectations
    for (int n = 0; n < 5; ++n) {
        EXPECT_NEAR(src_x_vals[n], coeff_golds[n] * golds[0], tol);
        EXPECT_NEAR(src_y_vals[n], coeff_golds[n] * golds[1], tol);
        EXPECT_NEAR(src_z_vals[n], coeff_golds[n] * golds[2], tol);
    }
}

TEST_F(ABLOffshoreMeshTest, boussinesq)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    // Initialize parameters
    populate_parameters();
    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    pde_mgr.register_transport_pde("Temperature");
    sim().create_transport_model();
    sim().init_physics();
    auto& mphase = sim().physics_manager().get<amr_wind::MultiPhase>();
    // Make sure to read water level
    mphase.post_init_actions();

    auto& src_term = pde_mgr.icns().fields().src_term;

    auto& temperature =
        sim().repo().get_field("temperature", amr_wind::FieldState::Old);
    auto& volume_fraction = sim().repo().get_field("vof");

    src_term.setVal(0.0_rt);
    auto& geom = sim().mesh().Geom();
    auto waterlev =
        sim().physics_manager().get<amr_wind::MultiPhase>().water_level();
    // Ensure temperature gradient crosses interface for the sake of testing
    const amrex::Real btm_temp_ht = -20;
    run_algorithm(temperature, [&](const int lev, const amrex::MFIter& mfi) {
        const auto bx = mfi.validbox();
        const auto& temp_arr = temperature(lev).array(mfi);
        init_abl_temperature_field(geom[lev], bx, temp_arr, btm_temp_ht);
        const auto& vof_arr = volume_fraction(lev).array(mfi);
        init_vof_field(geom[lev], bx, vof_arr, waterlev);
    });

    amr_wind::pde::icns::BoussinesqBuoyancy bb(sim());
    run_algorithm(temperature, [&](const int lev, const amrex::MFIter& mfi) {
        const auto bx = mfi.validbox();
        const auto& src_arr = src_term(lev).array(mfi);
        bb(lev, mfi, bx, amr_wind::FieldState::Old, src_arr);
    });

    // should be no forcing in x and y directions
    for (int i = 0; i < 2; ++i) {
        const auto min_src = utils::field_min(src_term, i);
        const auto max_src = utils::field_max(src_term, i);
        EXPECT_NEAR(min_src, 0.0_rt, tol);
        EXPECT_NEAR(max_src, 0.0_rt, tol);
    }

    // f = beta * (T0 - T)*g
    EXPECT_NEAR(utils::field_min(src_term, 2), 0.0_rt, tol);
    EXPECT_NEAR(
        utils::field_max(src_term, 2),
        -9.81_rt * (300.0_rt - 308.0_rt) / 300.0_rt, tol);

    // Check that the value is 0 below the interface but above btm_temp_ht
    const amrex::Real dz = geom[0].CellSize(2);
    const amrex::Real ploz = geom[0].ProbLo(2);
    // Divide by nx*ny cells because of sum
    const int nx = 8;
    const int ny = 8;
    amrex::Real src_z_val =
        get_val_at_height(
            src_term, 0, 0, ploz, dz, 0.5_rt * (waterlev + btm_temp_ht)) /
        nx / ny;
    // Check src value against expectations
    EXPECT_NEAR(src_z_val, 0.0_rt, tol);
}

} // namespace amr_wind_tests
