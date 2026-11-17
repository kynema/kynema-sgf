#include "gtest/gtest.h"
#include "src/physics/DampingLayer.H"
#include "src/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

TEST(DampingLayerMath, blending_function_type_parse)
{
    using namespace kynema_sgf::damping_layer;

    EXPECT_EQ(
        string_to_blending_function_type("linear"),
        BlendingFunctionType::Linear);
    EXPECT_EQ(
        string_to_blending_function_type("QUADRATIC"),
        BlendingFunctionType::Quadratic);
    EXPECT_EQ(
        string_to_blending_function_type("Exponential"),
        BlendingFunctionType::Exponential);
    EXPECT_EQ(
        string_to_blending_function_type("cosine"),
        BlendingFunctionType::Cosine);
}

TEST(DampingLayerMath, blending_function_values)
{
    using namespace kynema_sgf::damping_layer;

    constexpr amrex::Real tol = kynema_sgf::constants::TIGHT_TOL;

    EXPECT_NEAR(
        blending_function(0.0_rt, BlendingFunctionType::Linear), 1.0_rt, tol);
    EXPECT_NEAR(
        blending_function(1.0_rt, BlendingFunctionType::Linear), 0.0_rt, tol);

    EXPECT_NEAR(
        blending_function(0.5_rt, BlendingFunctionType::Quadratic), 0.75_rt,
        tol);

    const amrex::Real exp_mid =
        blending_function(0.5_rt, BlendingFunctionType::Exponential);
    EXPECT_GT(exp_mid, 0.0_rt);
    EXPECT_LT(exp_mid, 1.0_rt);

    EXPECT_NEAR(
        blending_function(0.0_rt, BlendingFunctionType::Cosine), 1.0_rt, tol);
    EXPECT_NEAR(
        blending_function(0.5_rt, BlendingFunctionType::Cosine), 0.5_rt, tol);
    EXPECT_NEAR(
        blending_function(1.0_rt, BlendingFunctionType::Cosine), 0.0_rt, tol);
}

TEST(DampingLayerMath, damping_calc_piecewise_behavior)
{
    using namespace kynema_sgf::damping_layer;

    constexpr amrex::Real tol = kynema_sgf::constants::TIGHT_TOL;
    constexpr amrex::Real thickness = 10.0_rt;
    constexpr amrex::Real blend_frac = 0.3_rt;
    constexpr amrex::Real full_damp_len = thickness * (1.0_rt - blend_frac);

    EXPECT_NEAR(
        damping_calc(
            2.0_rt, thickness, blend_frac, BlendingFunctionType::Cosine),
        1.0_rt, tol);

    const amrex::Real blend_pos =
        (8.0_rt - full_damp_len) / (thickness * blend_frac);
    const amrex::Real expected_mid =
        blending_function(blend_pos, BlendingFunctionType::Cosine);
    EXPECT_NEAR(
        damping_calc(
            8.0_rt, thickness, blend_frac, BlendingFunctionType::Cosine),
        expected_mid, tol);

    EXPECT_NEAR(
        damping_calc(
            11.0_rt, thickness, blend_frac, BlendingFunctionType::Cosine),
        0.0_rt, tol);

    EXPECT_NEAR(
        damping_calc(5.0_rt, thickness, 0.0_rt, BlendingFunctionType::Linear),
        1.0_rt, tol);
    EXPECT_NEAR(
        damping_calc(10.01_rt, thickness, 0.0_rt, BlendingFunctionType::Linear),
        0.0_rt, tol);
}

} // namespace kynema_sgf_tests
