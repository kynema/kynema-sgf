#include <AMReX_Gpu.H>
#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/eos_models/EOSModel.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

void test_eos_impl()
{
    const auto eos = amr_wind::eos::GammaLaw(1.01325e5_rt);
    amrex::Gpu::DeviceScalar<amrex::Real> val(0.0_rt);
    auto* d_val = val.dataPtr();

    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*unused*/) {
        d_val[0] = eos.p_rth(1.225_rt, 300.0_rt, 0.5_rt);
    });
    EXPECT_NEAR(
        val.dataValue(), 244859.65251771925_rt, amr_wind::constants::LOOSE_TOL);

    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*unused*/) {
        d_val[0] = eos.dp_constanttheta(1.225_rt, 300.0_rt, 0.5_rt);
    });
    EXPECT_NEAR(
        val.dataValue(), 279839.60287739342_rt, amr_wind::constants::LOOSE_TOL);
}

TEST(EOSModelTest, test_eos) { test_eos_impl(); }
} // namespace amr_wind_tests
