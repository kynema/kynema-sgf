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
    const amrex::Vector<amrex::Real> start{2.0_rt, 3.0_rt, 5.0_rt};
    amrex::Vector<amrex::Real> end{10.0_rt, 3.0_rt, 5.0_rt};

    // Test translation
    auto translate = amr_wind::channelbuilder::transform_to_local_coordinates(
        0.0_rt, 0.0_rt, 0.0_rt, start[0], start[1], start[2], end[0], end[1],
        end[2]);
    EXPECT_NEAR(translate[0], -2.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(translate[1], -3.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(translate[2], -5.0_rt, 1.0e-8_rt);

    // Test rotation
    end[0] = 7.0_rt;
    end[1] = 8.0_rt;
    end[2] = 10.0_rt;

    auto rotate = amr_wind::channelbuilder::transform_to_local_coordinates(
        end[0], end[1], end[2], start[0], start[1], start[2], end[0], end[1],
        end[2]);
    EXPECT_NEAR(rotate[0], std::sqrt(75.0_rt), 1.0e-8_rt);
    EXPECT_NEAR(rotate[1], 0.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(rotate[2], 0.0_rt, 1.0e-8_rt);

    // Test both in 2D - xy
    end[2] = 5.0_rt;

    auto both_xy = amr_wind::channelbuilder::transform_to_local_coordinates(
        2.0_rt + 2.0_rt / sqrt(2.0_rt), 3.0_rt, 5.0_rt, start[0], start[1], start[2], end[0], end[1],
        end[2]);
    EXPECT_NEAR(both_xy[0], 1.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(both_xy[1], -1.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(both_xy[2], 0.0_rt, 1.0e-8_rt);

    // Test both in 2D - xz
    end[1] = 3.0_rt;
    end[2] = 10.0_rt;

    auto both_xz = amr_wind::channelbuilder::transform_to_local_coordinates(
        2.0_rt + 2.0_rt / sqrt(2.0_rt), 3.0_rt, 5.0_rt, start[0], start[1], start[2], end[0], end[1],
        end[2]);
    EXPECT_NEAR(both_xz[0], 1.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(both_xz[1], 0.0_rt, 1.0e-8_rt);
    EXPECT_NEAR(both_xz[2], -1.0_rt, 1.0e-8_rt);

}

} // namespace amr_wind_tests
