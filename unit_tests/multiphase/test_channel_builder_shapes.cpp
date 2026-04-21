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

TEST(ChannelBuilderShapes, is_point_within_planes)
{
    // Define a simple segment from (0, 0, 0) to (10, 0, 0)
    const amrex::Real start_x = 0.0_rt;
    const amrex::Real start_y = 0.0_rt;
    const amrex::Real start_z = 0.0_rt;
    const amrex::Real end_x = 10.0_rt;
    const amrex::Real end_y = 0.0_rt;
    const amrex::Real end_z = 0.0_rt;

    // Point between planes should return true
    EXPECT_TRUE(amr_wind::channelbuilder::is_point_within_planes(
        5.0_rt, 0.0_rt, 0.0_rt, start_x, start_y, start_z, end_x, end_y,
        end_z));

    // Points on or near segment should return true
    EXPECT_TRUE(amr_wind::channelbuilder::is_point_within_planes(
        0.0_rt, 0.0_rt, 0.0_rt, start_x, start_y, start_z, end_x, end_y,
        end_z));
    EXPECT_TRUE(amr_wind::channelbuilder::is_point_within_planes(
        10.0_rt, 0.0_rt, 0.0_rt, start_x, start_y, start_z, end_x, end_y,
        end_z));

    // Points outside segment planes should return false
    EXPECT_FALSE(amr_wind::channelbuilder::is_point_within_planes(
        -5.0_rt, 0.0_rt, 0.0_rt, start_x, start_y, start_z, end_x, end_y,
        end_z));
    EXPECT_FALSE(amr_wind::channelbuilder::is_point_within_planes(
        15.0_rt, 0.0_rt, 0.0_rt, start_x, start_y, start_z, end_x, end_y,
        end_z));
}

TEST(ChannelBuilderShapes, transform_to_local_coordinates)
{
    // Define a simple segment from (0, 0, 0) to (10, 0, 0)
    // (along x-axis)
    const amrex::Real start_x = 0.0_rt;
    const amrex::Real start_y = 0.0_rt;
    const amrex::Real start_z = 0.0_rt;
    const amrex::Real end_x = 10.0_rt;
    const amrex::Real end_y = 0.0_rt;
    const amrex::Real end_z = 0.0_rt;

    // Point at segment start should have local coords near (0, 0, 0)
    auto local = amr_wind::channelbuilder::transform_to_local_coordinates(
        0.0_rt, 0.0_rt, 0.0_rt, start_x, start_y, start_z, end_x, end_y,
        end_z);
    EXPECT_NEAR(local[0], 0.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(local[1], 0.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(local[2], 0.0_rt, 1.0e-8_rt);

    // Point at segment midpoint should have xloc near 5
    local = amr_wind::channelbuilder::transform_to_local_coordinates(
        5.0_rt, 0.0_rt, 0.0_rt, start_x, start_y, start_z, end_x, end_y,
        end_z);
    EXPECT_NEAR(local[0], 5.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(local[1], 0.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(local[2], 0.0_rt, 1.0e-8_rt);

    // Point offset in y should have yloc with offset value
    local = amr_wind::channelbuilder::transform_to_local_coordinates(
        5.0_rt, 2.0_rt, 0.0_rt, start_x, start_y, start_z, end_x, end_y,
        end_z);
    EXPECT_NEAR(local[0], 5.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(local[1], 2.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(local[2], 0.0_rt, 1.0e-8_rt);

    // Point offset in z should have zloc with offset value
    local = amr_wind::channelbuilder::transform_to_local_coordinates(
        5.0_rt, 0.0_rt, 3.0_rt, start_x, start_y, start_z, end_x, end_y,
        end_z);
    EXPECT_NEAR(local[0], 5.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(local[1], 0.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(local[2], 3.0_rt, 1.0e-8_rt);
}

TEST(ChannelBuilderShapes, transform_to_local_coordinates_diagonal)
{
    // Define a diagonal segment from (0, 0, 0) to (10, 10, 0)
    // This tests rotation in the xy plane
    const amrex::Real start_x = 0.0_rt;
    const amrex::Real start_y = 0.0_rt;
    const amrex::Real start_z = 0.0_rt;
    const amrex::Real end_x = 10.0_rt;
    const amrex::Real end_y = 10.0_rt;
    const amrex::Real end_z = 0.0_rt;

    // Point at segment start
    auto local = amr_wind::channelbuilder::transform_to_local_coordinates(
        0.0_rt, 0.0_rt, 0.0_rt, start_x, start_y, start_z, end_x, end_y,
        end_z);
    EXPECT_NEAR(local[0], 0.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(local[1], 0.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(local[2], 0.0_rt, 1.0e-8_rt);

    // Point at segment in the direction of (end - start)
    // Should have xloc > 0 and yloc ~= 0
    local = amr_wind::channelbuilder::transform_to_local_coordinates(
        7.071067812_rt, 7.071067812_rt, 0.0_rt, start_x, start_y, start_z,
        end_x, end_y, end_z);
    EXPECT_NEAR(local[0], 10.0_rt, 1.0e-5_rt);
    EXPECT_NEAR(local[1], 0.0_rt, 1.0e-5_rt);
    EXPECT_NEAR(local[2], 0.0_rt, 1.0e-5_rt);
}

} // namespace amr_wind_tests
