#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/core/vs/vector_space.H"
#include "amr-wind/wind_energy/actuator/actuator_utils.H"
#include <cmath>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace act = ::amr_wind::actuator::utils;
namespace vs = ::amr_wind::vs;
namespace utils = ::amr_wind::utils;

namespace amr_wind_tests::amr_wind {
namespace {

TEST(CylindricalTransformation, distances_on_z_aligned_circle)
{
    const vs::Vector p1{0, 2, 1};
    const vs::Vector p2{1, 0, 0};
    const vs::Vector normal{0, 0, 1};
    const vs::Vector origin{0, 0, 0};

    auto d = act::delta_pnts_cyl(origin, normal, p1, p2);
    EXPECT_DOUBLE_EQ(1.0_rt, d[0]);
    EXPECT_DOUBLE_EQ(utils::half_pi(), d[1]);
    EXPECT_DOUBLE_EQ(1.0_rt, d[2]);
}

TEST(CylindricalTransformation, returns_correct_radius_when_origin_is_a_point)
{
    const vs::Vector p2{0, 1, 0};
    const vs::Vector normal{0, 0, 1};
    const vs::Vector origin{0, 0, 0};
    const vs::Vector p1(origin);

    auto d = act::delta_pnts_cyl(origin, normal, p1, p2);
    EXPECT_DOUBLE_EQ(1.0_rt, d[0]);
    // theta will be ill defined
    EXPECT_DOUBLE_EQ(0.0_rt, d[2]);
}

TEST(CylindricalTransformation, distances_on_shifted_circle)
{
    const vs::Vector p1{1, 3, 2};
    const vs::Vector p2{2, 1, 1};
    const vs::Vector normal{0, 0, 1};
    const vs::Vector origin{1, 1, 1};

    auto d = act::delta_pnts_cyl(origin, normal, p1, p2);
    EXPECT_DOUBLE_EQ(1.0_rt, d[0]);
    EXPECT_DOUBLE_EQ(utils::half_pi(), d[1]);
    EXPECT_DOUBLE_EQ(1.0_rt, d[2]);
}

vs::Vector rotation(const vs::Vector& angles, const vs::Vector& data)
{
    const vs::Tensor rotMatrix =
        vs::xrot(angles.x()) & vs::yrot(angles.y()) & vs::zrot(angles.z());
    return data & rotMatrix;
}

TEST(CylindricalTransformation, distances_on_rotated_circle)
{
    const vs::Vector p1{0, 2, 1};
    const vs::Vector p2{1, 0, 0};
    const vs::Vector normal{0, 0, 1};
    const vs::Vector origin{0, 0, 0};
    const vs::Vector angles{45, -80, 10};

    auto newNormal = rotation(angles, normal);
    auto np1 = rotation(angles, p1);
    auto np2 = rotation(angles, p2);

    auto d = act::delta_pnts_cyl(origin, newNormal, np1, np2);
    EXPECT_DOUBLE_EQ(1.0_rt, d[0]);
    EXPECT_DOUBLE_EQ(utils::half_pi(), d[1])
        << "angle in degress: " << utils::degrees(d[1]);
    EXPECT_NEAR(1.0_rt, d[2], std::numeric_limits<float>::epsilon() * 1.0e1_rt);
}

TEST(CylindricalTransformation, distances_on_rotated_shifted_circle)
{
    // for this test we define the points relative to the origin to make
    // it easier to use an analytic solution
    const vs::Vector p1{0, 2, 1};
    const vs::Vector p2{1, 0, 0};
    const vs::Vector normal{0, 0, 1};
    const vs::Vector origin{5, 5, 5};
    const vs::Vector angles{45, -80, 10};

    // perform rotations about the Cartesian coordinate system origin to define
    // a new normal and points relative to that normal then add the origin to
    // shift the point locations relative to the origin
    auto newNormal = rotation(angles, normal);
    auto np1 = rotation(angles, p1) + origin;
    auto np2 = rotation(angles, p2) + origin;

    auto d = act::delta_pnts_cyl(origin, newNormal, np1, np2);
    EXPECT_NEAR(
        1.0_rt, d[0], std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
    EXPECT_NEAR(
        utils::half_pi(), d[1],
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt)
        << "angle in degress: " << utils::degrees(d[1]);
    EXPECT_NEAR(
        1.0_rt, d[2], std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
}

} // namespace
} // namespace amr_wind_tests::amr_wind
