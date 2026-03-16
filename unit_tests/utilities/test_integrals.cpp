#include <AMReX_Gpu.H>
#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/utilities/integrals.H"
#include "amr-wind/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

void test_trapezoid_integration_xsquared_impl()
{
    const amrex::Real xa = 0.0_rt;
    const amrex::Real xb = 1.2_rt;
    const int n = 10000;

    amrex::Gpu::DeviceScalar<amrex::Real> integ(0.0_rt);
    auto* d_integ = integ.dataPtr();
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*unused*/) {
        d_integ[0] = amr_wind::utils::trapz(
            xa, xb, n, [](const amrex::Real x) { return x * x; });
    });

    EXPECT_NEAR(integ.dataValue(), 0.576_rt, amr_wind::constants::LOOSE_TOL);
}

TEST(Integrals, trapezoid_integration_xsquared)
{
    test_trapezoid_integration_xsquared_impl();
}

} // namespace amr_wind_tests
