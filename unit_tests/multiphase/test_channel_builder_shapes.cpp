#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/physics/multiphase/ChannelBuilder.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

TEST(ChannelBuilderShapes, trapezoid_inside_outside)
{
    const amrex::Real top = 2.0_rt;
    const amrex::Real bottom = 4.0_rt;
    const amrex::Real height = 6.0_rt;

    // Interior points
    EXPECT_TRUE(amr_wind::channelbuilder::trapezoid(
        top, bottom, height, 0.0_rt, 0.0_rt));
    EXPECT_TRUE(amr_wind::channelbuilder::trapezoid(
        top, bottom, height, 0.8_rt, 1.5_rt));

    // Boundary point should be included
    EXPECT_TRUE(amr_wind::channelbuilder::trapezoid(
        top, bottom, height, -1.5_rt, 3.0_rt));

    // Exterior points
    EXPECT_FALSE(amr_wind::channelbuilder::trapezoid(
        top, bottom, height, 3.1_rt, 0.0_rt));
    EXPECT_FALSE(amr_wind::channelbuilder::trapezoid(
        top, bottom, height, 0.0_rt, 3.1_rt));
}

TEST(ChannelBuilderShapes, ellipse_inside_outside)
{
    const amrex::Real ax_horz = 3.0_rt;
    const amrex::Real ax_vert = 2.0_rt;

    // Interior points
    EXPECT_TRUE(amr_wind::channelbuilder::ellipse(
        ax_horz, ax_vert, 0.0_rt, 0.0_rt));
    EXPECT_TRUE(amr_wind::channelbuilder::ellipse(
        ax_horz, ax_vert, 2.0_rt, 0.5_rt));

    // Boundary point should be included
    EXPECT_TRUE(amr_wind::channelbuilder::ellipse(
        ax_horz, ax_vert, 3.0_rt, 0.0_rt));

    // Exterior points
    EXPECT_FALSE(amr_wind::channelbuilder::ellipse(
        ax_horz, ax_vert, 3.2_rt, 0.0_rt));
    EXPECT_FALSE(amr_wind::channelbuilder::ellipse(
        ax_horz, ax_vert, 2.5_rt, 1.5_rt));
}

} // namespace amr_wind_tests
