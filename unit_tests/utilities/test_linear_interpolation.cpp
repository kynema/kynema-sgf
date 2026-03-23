#include <algorithm>
#include <numbers>
#include <utility>
#include <vector>
#include <numeric>
#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_Random.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

TEST(LinearInterpolation, check_bounds)
{
    namespace interp = amr_wind::interp;
    std::vector<amrex::Real> xvec(10);
    std::iota(xvec.begin(), xvec.end(), 0.0_rt);

    const auto* start = xvec.data();
    const auto* end = xvec.data() + xvec.size();
    {
        const auto idx = interp::check_bounds(start, end, -1.0_rt);
        EXPECT_EQ(idx.idx, 0);
        EXPECT_EQ(idx.lim, interp::Limits::LOWLIM);
    }
    {
        const auto idx = interp::check_bounds(start, end, 9.1_rt);
        EXPECT_EQ(idx.idx, 9);
        EXPECT_EQ(idx.lim, interp::Limits::UPLIM);
    }
    {
        const amrex::Real xinp = 9.0_rt * amrex::Random();
        const auto idx = interp::check_bounds(start, end, xinp);
        EXPECT_EQ(idx.idx, 0);
        EXPECT_EQ(idx.lim, interp::Limits::VALID);
    }
}

TEST(LinearInterpolation, bisection_search)
{
    namespace interp = amr_wind::interp;
    std::vector<amrex::Real> xvec(10);
    std::iota(xvec.begin(), xvec.end(), 0.0_rt);

    const auto* start = xvec.data();
    const auto* end = xvec.data() + xvec.size();
    {
        const auto idx = interp::bisection_search(start, end, 5.0_rt);
        EXPECT_EQ(idx.idx, 4);
        EXPECT_EQ(idx.lim, interp::Limits::VALID);
    }
    {
        const auto idx = interp::bisection_search(start, end, -1.0_rt);
        EXPECT_EQ(idx.idx, 0);
        EXPECT_EQ(idx.lim, interp::Limits::LOWLIM);
    }
    {
        const auto idx = interp::bisection_search(start, end, 9.1_rt);
        EXPECT_EQ(idx.idx, 9);
        EXPECT_EQ(idx.lim, interp::Limits::UPLIM);
    }
}

TEST(LinearInterpolation, nearest_search)
{
    namespace interp = amr_wind::interp;
    std::vector<amrex::Real> xvec(10);
    std::iota(xvec.begin(), xvec.end(), 0.0_rt);

    const auto* start = xvec.data();
    const auto* end = xvec.data() + xvec.size();
    {
        const auto idx = interp::nearest_search(start, end, 5.0_rt);
        EXPECT_EQ(idx.idx, 5);
        EXPECT_EQ(idx.lim, interp::Limits::VALID);
    }
    {
        const auto idx = interp::nearest_search(start, end, 1.1_rt);
        EXPECT_EQ(idx.idx, 1);
        EXPECT_EQ(idx.lim, interp::Limits::VALID);
    }
    {
        const auto idx = interp::nearest_search(start, end, 1.6_rt);
        EXPECT_EQ(idx.idx, 2);
        EXPECT_EQ(idx.lim, interp::Limits::VALID);
    }
    {
        const auto idx = interp::nearest_search(start, end, 3.5_rt);
        EXPECT_EQ(idx.idx, 3);
        EXPECT_EQ(idx.lim, interp::Limits::VALID);
    }
    {
        const auto idx = interp::nearest_search(start, end, -1.0_rt);
        EXPECT_EQ(idx.idx, 0);
        EXPECT_EQ(idx.lim, interp::Limits::LOWLIM);
    }
    {
        const auto idx = interp::nearest_search(start, end, 9.1_rt);
        EXPECT_EQ(idx.idx, 9);
        EXPECT_EQ(idx.lim, interp::Limits::UPLIM);
    }
}

TEST(LinearInterpolation, find_index)
{
    namespace interp = amr_wind::interp;
    std::vector<amrex::Real> xvec(10);
    std::iota(xvec.begin(), xvec.end(), 0.0_rt);

    const auto* start = xvec.data();
    const auto* end = xvec.data() + xvec.size();
    {
        const auto idx = interp::find_index(start, end, 5.0_rt);
        EXPECT_EQ(idx.idx, 4);
        EXPECT_EQ(idx.lim, interp::Limits::VALID);
    }
    {
        const auto idx = interp::find_index(start, end, 5.0_rt, 3);
        EXPECT_EQ(idx.idx, 4);
        EXPECT_EQ(idx.lim, interp::Limits::VALID);
    }
    {
        const auto idx = interp::find_index(start, end, -1.0_rt);
        EXPECT_EQ(idx.idx, 0);
        EXPECT_EQ(idx.lim, interp::Limits::LOWLIM);
    }
    {
        const auto idx = interp::find_index(start, end, 9.1_rt);
        EXPECT_EQ(idx.idx, 9);
        EXPECT_EQ(idx.lim, interp::Limits::UPLIM);
    }
}

TEST(LinearInterpolation, lin_interp_single)
{
    namespace interp = amr_wind::interp;

    const amrex::Real mult_fac = 2.0_rt + (10.0_rt * amrex::Random());
    std::vector<amrex::Real> xvec(10), yvec(10);
    std::iota(xvec.begin(), xvec.end(), 0.0_rt);
    std::ranges::transform(xvec, yvec.begin(), [mult_fac](const amrex::Real x) {
        return mult_fac * x;
    });

    std::vector<amrex::Real> xtest{2.5_rt, 4.5_rt, 6.3_rt, 8.8_rt};
    for (const auto& x : xtest) {
        const auto y = interp::linear(xvec, yvec, x);
        EXPECT_NEAR(
            y, mult_fac * x,
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
    }
}

TEST(LinearInterpolation, lin_interp_single_multicomponent)
{
    namespace interp = amr_wind::interp;

    const int ncomp = 3;
    std::vector<amrex::Real> xvec(10);
    std::vector<amrex::Real> yvec(static_cast<unsigned long>(ncomp * 10));
    std::iota(xvec.begin(), xvec.end(), 0.0_rt);
    const amrex::Vector<amrex::Real> mult_facs = {
        2.0_rt + (10.0_rt * amrex::Random()),
        2.0_rt + (10.0_rt * amrex::Random()),
        2.0_rt + (10.0_rt * amrex::Random())};
    for (int i = 0; std::cmp_less(i, xvec.size()); i++) {
        for (int n = 0; n < ncomp; n++) {
            yvec[(ncomp * i) + n] = mult_facs[n] * xvec[i];
        }
    }

    std::vector<amrex::Real> xtest{2.5_rt, 4.5_rt, 6.3_rt, 8.8_rt};
    for (const auto& x : xtest) {
        for (int n = 0; n < ncomp; n++) {
            const auto y = interp::linear(xvec, yvec, x, 3, n);
            EXPECT_NEAR(
                y, mult_facs[n] * x,
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
        }
    }
}

TEST(LinearInterpolation, bilin_interp_single)
{
    namespace interp = amr_wind::interp;

    const amrex::Real mult_facx = 2.0_rt + (10.0_rt * amrex::Random());
    const amrex::Real mult_facy = 2.0_rt + (10.0_rt * amrex::Random());
    std::vector<amrex::Real> xvec(10), yvec(5);
    std::vector<amrex::Real> zvec(xvec.size() * yvec.size());
    std::iota(xvec.begin(), xvec.end(), 0.0_rt);
    std::iota(yvec.begin(), yvec.end(), 0.0_rt);

    for (int i = 0; std::cmp_less(i, xvec.size()); i++) {
        for (int j = 0; std::cmp_less(j, yvec.size()); j++) {
            zvec[(i * yvec.size()) + j] =
                mult_facx * xvec[i] * mult_facy * yvec[j];
        }
    }

    std::vector<amrex::Real> xtest{2.5_rt, 4.5_rt, 6.3_rt, 8.8_rt};
    std::vector<amrex::Real> ytest{1.1_rt, 2.3_rt, 3.1_rt, 3.8_rt};
    for (const auto& x : xtest) {
        for (const auto& y : ytest) {
            const auto z = interp::bilinear(xvec, yvec, zvec, x, y);
            EXPECT_NEAR(
                z, mult_facx * x * mult_facy * y,
                std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
        }
    }
}

TEST(LinearInterpolation, lin_interp)
{
    namespace interp = amr_wind::interp;

    const amrex::Real mult_fac = 2.0_rt + (10.0_rt * amrex::Random());
    std::vector<amrex::Real> xvec(10), yvec(10);
    std::iota(xvec.begin(), xvec.end(), 0.0_rt);
    std::ranges::transform(xvec, yvec.begin(), [mult_fac](const amrex::Real x) {
        return mult_fac * x;
    });

    std::vector<amrex::Real> xtest{2.5_rt, 4.5_rt, 6.3_rt, 8.8_rt};
    std::vector<amrex::Real> ytest(xtest.size());

    interp::linear(xvec, yvec, xtest, ytest);
    for (size_t i = 0; i < xtest.size(); ++i) {
        EXPECT_NEAR(
            ytest[i], mult_fac * xtest[i],
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
    }
}

TEST(LinearInterpolation, lin_monotonic)
{
    namespace interp = amr_wind::interp;

    const amrex::Real mult_fac = 2.0_rt + (10.0_rt * amrex::Random());
    std::vector<amrex::Real> xvec(10), yvec(10);
    std::iota(xvec.begin(), xvec.end(), 0.0_rt);
    std::ranges::transform(xvec, yvec.begin(), [mult_fac](const amrex::Real x) {
        return mult_fac * x;
    });

    std::vector<amrex::Real> xtest{2.5_rt, 4.5_rt, 6.3_rt, 8.8_rt};
    std::vector<amrex::Real> ytest(xtest.size());

    interp::linear_monotonic(xvec, yvec, xtest, ytest);
    for (size_t i = 0; i < xtest.size(); ++i) {
        EXPECT_NEAR(
            ytest[i], mult_fac * xtest[i],
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
    }
}

TEST(LinearInterpolation, lin_interp_angle)
{
    namespace interp = amr_wind::interp;

    const int vecsize = 8;
    std::vector<amrex::Real> xvec(vecsize);
    std::iota(xvec.begin(), xvec.end(), 0.0_rt);
    // Set up with angles that are awkward to interpolate between
    std::vector<amrex::Real> yvec_deg{0.0_rt,   30.0_rt,   330.0_rt, 300.0_rt,
                                      -15.0_rt, -180.0_rt, 120.0_rt, -150.0_rt};
    // Create duplicate list in radians
    std::vector<amrex::Real> yvec_rad(vecsize);
    for (size_t i = 0; i < yvec_rad.size(); ++i) {
        yvec_rad[i] = amr_wind::utils::radians(yvec_deg[i]);
    }

    // Set up output vectors, interp locations, and golds
    std::vector<amrex::Real> xtest(vecsize), ytest_deg(vecsize),
        ytest_rad(vecsize);
    std::iota(xtest.begin(), xtest.end(), 0.5_rt);
    xtest[7] = 4.0_rt;
    std::vector<amrex::Real> ygold_deg{
        15.0_rt,  0.0_rt,    315.0_rt, 322.5_rt,
        -97.5_rt, -210.0_rt, 165.0_rt, 360.0_rt - 15.0_rt};

    // Upper bound is 360 because y is in degrees
    interp::linear_angle(xvec, yvec_deg, xtest, ytest_deg, 360.0_rt);
    // Upper bound is 2pi because y is in radians
    interp::linear_angle(
        xvec, yvec_rad, xtest, ytest_rad,
        2.0_rt * std::numbers::pi_v<amrex::Real>);
    for (size_t i = 0; i < xtest.size(); ++i) {
        EXPECT_NEAR(
            ytest_deg[i], ygold_deg[i],
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
        EXPECT_NEAR(
            ytest_rad[i], amr_wind::utils::radians(ygold_deg[i]),
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
    }
}

} // namespace amr_wind_tests
