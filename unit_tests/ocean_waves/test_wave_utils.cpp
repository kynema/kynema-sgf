#include <algorithm>
#include <numbers>
#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/ocean_waves/utils/wave_utils_K.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

class WaveUtilsTest : public MeshTest
{};

TEST_F(WaveUtilsTest, free_surface_to_vof)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    constexpr amrex::Real zsl = 1.0_rt;
    constexpr amrex::Real dz = 0.5_rt;

    const amrex::Real vof_fully_above =
        amr_wind::ocean_waves::utils::free_surface_to_vof(zsl, 100.0_rt, dz);
    const amrex::Real vof_fully_below =
        amr_wind::ocean_waves::utils::free_surface_to_vof(zsl, -100.0_rt, dz);
    const amrex::Real vof_at_level =
        amr_wind::ocean_waves::utils::free_surface_to_vof(zsl, zsl, dz);
    const amrex::Real vof_half_above =
        amr_wind::ocean_waves::utils::free_surface_to_vof(
            zsl, zsl + (0.25_rt * dz), dz);
    const amrex::Real vof_half_below =
        amr_wind::ocean_waves::utils::free_surface_to_vof(
            zsl, zsl - (0.25_rt * dz), dz);
    const amrex::Real vof_exactly_above =
        amr_wind::ocean_waves::utils::free_surface_to_vof(
            zsl, zsl + (0.5_rt * dz), dz);
    const amrex::Real vof_exactly_below =
        amr_wind::ocean_waves::utils::free_surface_to_vof(
            zsl, zsl - (0.5_rt * dz), dz);

    EXPECT_NEAR(vof_fully_above, 0.0_rt, tol);
    EXPECT_NEAR(vof_fully_below, 1.0_rt, tol);
    EXPECT_NEAR(vof_at_level, 0.5_rt, tol);
    EXPECT_NEAR(vof_half_above, 0.25_rt, tol);
    EXPECT_NEAR(vof_half_below, 0.75_rt, tol);
    EXPECT_NEAR(vof_exactly_above, 0.0_rt, tol);
    EXPECT_NEAR(vof_exactly_below, 1.0_rt, tol);
}

TEST_F(WaveUtilsTest, gamma_generate)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    constexpr amrex::Real zone_length = 1.5_rt;

    const amrex::Real gamma_past_left =
        amr_wind::ocean_waves::utils::gamma_generate(-1.0_rt, zone_length);
    const amrex::Real gamma_left =
        amr_wind::ocean_waves::utils::gamma_generate(0.0_rt, zone_length);
    const amrex::Real gamma_middle =
        amr_wind::ocean_waves::utils::gamma_generate(
            0.5_rt * zone_length, zone_length);
    const amrex::Real gamma_right =
        amr_wind::ocean_waves::utils::gamma_generate(zone_length, zone_length);
    const amrex::Real gamma_past_right =
        amr_wind::ocean_waves::utils::gamma_generate(
            zone_length + 1.0_rt, zone_length);

    const amrex::Real gamma_middle_gold =
        1.0_rt - ((std::exp(std::pow(0.5_rt, 3.5_rt)) - 1.0_rt) /
                  (std::exp(1.0_rt) - 1.0_rt));

    EXPECT_NEAR(gamma_past_left, 0.0_rt, tol);
    EXPECT_NEAR(gamma_left, 0.0_rt, tol);
    EXPECT_NEAR(gamma_middle, gamma_middle_gold, tol);
    EXPECT_NEAR(gamma_right, 1.0_rt, tol);
    EXPECT_NEAR(gamma_past_right, 1.0_rt, tol);
}

TEST_F(WaveUtilsTest, gamma_absorb)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    constexpr amrex::Real zone_length = 1.5_rt;

    const amrex::Real gamma_past_left =
        amr_wind::ocean_waves::utils::gamma_absorb(
            -1.0_rt, zone_length, 1.0_rt);
    const amrex::Real gamma_left =
        amr_wind::ocean_waves::utils::gamma_absorb(0.0_rt, zone_length, 1.0_rt);
    const amrex::Real gamma_middle = amr_wind::ocean_waves::utils::gamma_absorb(
        0.5_rt * zone_length, zone_length, 1.0_rt);
    const amrex::Real gamma_right = amr_wind::ocean_waves::utils::gamma_absorb(
        zone_length, zone_length, 1.0_rt);
    const amrex::Real gamma_past_right =
        amr_wind::ocean_waves::utils::gamma_absorb(
            zone_length + 1.0_rt, zone_length, 1.0_rt);

    const amrex::Real gamma_right_length_factor2 =
        amr_wind::ocean_waves::utils::gamma_absorb(
            zone_length, zone_length, 2.0_rt);

    const amrex::Real gamma_middle_gold =
        1.0_rt - ((std::exp(std::pow(0.5_rt, 3.5_rt)) - 1.0_rt) /
                  (std::exp(1.0_rt) - 1.0_rt));

    EXPECT_NEAR(gamma_past_left, 1.0_rt, tol);
    EXPECT_NEAR(gamma_left, 1.0_rt, tol);
    EXPECT_NEAR(gamma_middle, gamma_middle_gold, tol);
    EXPECT_NEAR(gamma_right_length_factor2, gamma_middle_gold, tol);
    EXPECT_NEAR(gamma_right, 0.0_rt, tol);
    EXPECT_NEAR(gamma_past_right, 0.0_rt, tol);
}

TEST_F(WaveUtilsTest, ramp)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    constexpr amrex::Real near_tol = 1.0e-5_rt;
    constexpr amrex::Real ramp_period = 1.5_rt;

    const amrex::Real f_ramp_begin =
        amr_wind::ocean_waves::utils::ramp(0.0_rt, ramp_period);
    const amrex::Real f_ramp_middle =
        amr_wind::ocean_waves::utils::ramp(0.5_rt * ramp_period, ramp_period);
    const amrex::Real f_ramp_near_end =
        amr_wind::ocean_waves::utils::ramp(ramp_period - near_tol, ramp_period);
    const amrex::Real f_ramp_end =
        amr_wind::ocean_waves::utils::ramp(ramp_period, ramp_period);
    const amrex::Real f_ramp_past =
        amr_wind::ocean_waves::utils::ramp(ramp_period + 1.0_rt, ramp_period);

    const amrex::Real f_ramp_middle_gold =
        0.5_rt - (std::sin(std::numbers::pi_v<amrex::Real> * 0.5_rt) /
                  std::numbers::pi_v<amrex::Real>);

    EXPECT_NEAR(f_ramp_begin, 0.0_rt, tol);
    EXPECT_NEAR(f_ramp_middle, f_ramp_middle_gold, tol);
    EXPECT_NEAR(f_ramp_near_end, 1.0_rt, 2.0_rt * near_tol);
    EXPECT_NEAR(f_ramp_end, 1.0_rt, tol);
    EXPECT_NEAR(f_ramp_past, 1.0_rt, tol);
}

TEST_F(WaveUtilsTest, combine_linear)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    constexpr amrex::Real Gamma = 0.7_rt;
    constexpr amrex::Real target = 2.0_rt;
    constexpr amrex::Real current = 1.5_rt;

    const amrex::Real result =
        amr_wind::ocean_waves::utils::combine_linear(Gamma, target, current);

    EXPECT_NEAR(result, ((1.0_rt - Gamma) * target) + (Gamma * current), tol);
}

TEST_F(WaveUtilsTest, harmonize_profiles)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    constexpr amrex::Real problo_x = -1.0_rt;
    constexpr amrex::Real probhi_x = 1.0_rt;
    constexpr amrex::Real gen_length = 0.2_rt;
    constexpr amrex::Real beach_length = 0.4_rt;

    const amr_wind::ocean_waves::utils::WaveVec left{
        1.0_rt, 2.0_rt, 3.0_rt, 0.1_rt};
    const amr_wind::ocean_waves::utils::WaveVec bulk{
        0.0_rt, 1.1_rt, 2.0_rt, 0.0_rt};
    const amr_wind::ocean_waves::utils::WaveVec right{
        -1.0_rt, 1.0_rt, 0.5_rt, -0.1_rt};

    amrex::Real x = -0.9_rt;
    auto result = amr_wind::ocean_waves::utils::harmonize_profiles_1d(
        x, problo_x, gen_length, probhi_x, beach_length, left, bulk, right);
    for (int n = 0; n < 4; ++n) {
        EXPECT_NEAR(left[n], result[n], tol);
    }

    x = -0.75_rt;
    result = amr_wind::ocean_waves::utils::harmonize_profiles_1d(
        x, problo_x, gen_length, probhi_x, beach_length, left, bulk, right);
    const amrex::Real Gamma_l = amr_wind::ocean_waves::utils::gamma_generate(
        x - (problo_x + gen_length), 0.5_rt * gen_length);
    for (int n = 0; n < 4; ++n) {
        EXPECT_NEAR(
            ((1.0_rt - Gamma_l) * left[n]) + (Gamma_l * bulk[n]), result[n],
            tol);
    }

    x = 0.;
    result = amr_wind::ocean_waves::utils::harmonize_profiles_1d(
        x, problo_x, gen_length, probhi_x, beach_length, left, bulk, right);
    for (int n = 0; n < 4; ++n) {
        EXPECT_NEAR(bulk[n], result[n], tol);
    }

    x = 0.5_rt;
    result = amr_wind::ocean_waves::utils::harmonize_profiles_1d(
        x, problo_x, gen_length, probhi_x, beach_length, left, bulk, right);
    const amrex::Real Gamma_r = amr_wind::ocean_waves::utils::gamma_absorb(
        x - (probhi_x - beach_length) + (0.5_rt * beach_length),
        0.5_rt * beach_length, 1.0_rt);
    for (int n = 0; n < 4; ++n) {
        EXPECT_NEAR(
            ((1.0_rt - Gamma_r) * right[n]) + (Gamma_r * bulk[n]), result[n],
            tol);
    }

    x = 0.8_rt;
    result = amr_wind::ocean_waves::utils::harmonize_profiles_1d(
        x, problo_x, gen_length, probhi_x, beach_length, left, bulk, right);
    for (int n = 0; n < 4; ++n) {
        EXPECT_NEAR(right[n], result[n], tol);
    }
}

// Proof of concept test for multiple gammas in x and y
TEST_F(WaveUtilsTest, gamma_xy)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    constexpr amrex::Real zone_length = 0.25_rt;

    const amrex::Vector<amrex::Real> x{0.0_rt,
                                       1.0_rt,
                                       0.5_rt,
                                       0.5_rt,
                                       zone_length,
                                       zone_length,
                                       0.5_rt * zone_length};
    const amrex::Vector<amrex::Real> y{0.5_rt,
                                       0.5_rt,
                                       0.0_rt,
                                       1.0_rt,
                                       zone_length,
                                       0.5_rt,
                                       0.5_rt * zone_length};
    const amrex::Real gm =
        1.0_rt -
        ((std::exp(std::pow(1.0_rt - (x[x.size() - 1] / zone_length), 3.5_rt)) -
          1.0_rt) /
         (std::exp(1.0_rt) - 1.0_rt));
    const amrex::Vector<amrex::Real> gold_gamma{0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt,
                                                1.0_rt, 1.0_rt, gm};

    for (int n = 0; n < x.size(); ++n) {
        const amrex::Real gamma_xlo =
            amr_wind::ocean_waves::utils::gamma_generate(
                x[n] - 0.0_rt, zone_length);
        const amrex::Real gamma_xhi =
            amr_wind::ocean_waves::utils::gamma_absorb(
                x[n] - (1.0_rt - zone_length), zone_length, 1.0_rt);
        const amrex::Real gamma_ylo =
            amr_wind::ocean_waves::utils::gamma_generate(
                y[n] - 0.0_rt, zone_length);
        const amrex::Real gamma_yhi =
            amr_wind::ocean_waves::utils::gamma_absorb(
                y[n] - (1.0_rt - zone_length), zone_length, 1.0_rt);
        const amrex::Real gamma =
            std::min({gamma_xhi, gamma_xlo, gamma_yhi, gamma_ylo});
        EXPECT_NEAR(gamma, gold_gamma[n], tol);
    }
}

} // namespace amr_wind_tests
