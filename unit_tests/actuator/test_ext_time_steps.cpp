#include <string>
#include "ks_test_utils/AmrexTest.H"
#include "src/actuator/turbine/external/ExtTurbIface.H"
#include "AMReX_Exception.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {
namespace {

// Time steps of an external solver in a CFD time step
class ExtTimeStepsTest : public AmrexTest
{
protected:
    static int num_steps(const amrex::Real interval, const amrex::Real dt)
    {
        return ::ext_turb::num_whole_steps(
            interval, dt, "CFD time step", "solver time steps");
    }
};

} // namespace

TEST_F(ExtTimeStepsTest, ratio_just_below_whole_number)
{
    // Double: 0.3 / 0.1 = 2.9999999999999996, 0.35 / 0.05 = 6.999999999999999
    EXPECT_EQ(num_steps(0.3_rt, 0.1_rt), 3);
    EXPECT_EQ(num_steps(0.0003_rt, 0.0001_rt), 3);
    EXPECT_EQ(num_steps(0.009_rt, 0.003_rt), 3);
    EXPECT_EQ(num_steps(0.35_rt, 0.05_rt), 7);
}

TEST_F(ExtTimeStepsTest, ratio_just_above_whole_number)
{
    // Double: 0.035 / 0.005 = 7.000000000000001
    EXPECT_EQ(num_steps(0.035_rt, 0.005_rt), 7);
    EXPECT_EQ(num_steps(0.00875_rt, 0.00125_rt), 7);
}

TEST_F(ExtTimeStepsTest, whole_ratio)
{
    EXPECT_EQ(num_steps(0.00125_rt, 0.00125_rt), 1);
    EXPECT_EQ(num_steps(0.005_rt, 0.00125_rt), 4);
    EXPECT_EQ(num_steps(0.1_rt, 0.00125_rt), 80);
    // Replay from the initial time: no step
    EXPECT_EQ(num_steps(0.0_rt, 0.00125_rt), 0);
}

TEST_F(ExtTimeStepsTest, rejects_ratio_that_is_not_whole)
{
    EXPECT_THROW(num_steps(0.00125_rt, 0.001_rt), amrex::RuntimeError);
    // 3.0003: beyond the round-off tolerance in single and double
    EXPECT_THROW(num_steps(0.0030003_rt, 0.001_rt), amrex::RuntimeError);
    // Solver step larger than the CFD step
    EXPECT_THROW(num_steps(0.0004_rt, 0.001_rt), amrex::RuntimeError);
    EXPECT_THROW(num_steps(0.0006_rt, 0.001_rt), amrex::RuntimeError);
}

TEST_F(ExtTimeStepsTest, rejects_time_step_that_is_not_positive)
{
    EXPECT_THROW(num_steps(0.005_rt, 0.0_rt), amrex::RuntimeError);
    EXPECT_THROW(num_steps(0.005_rt, -0.00125_rt), amrex::RuntimeError);
    EXPECT_THROW(num_steps(-0.005_rt, 0.00125_rt), amrex::RuntimeError);
}

TEST_F(ExtTimeStepsTest, step_tolerance)
{
    const amrex::Real dt = 0.001_rt;
#ifndef AMREX_USE_FLOAT
    // Replay start time off by 2e-4 of a CFD step after 143000 steps
    const amrex::Real t_start = (143000.0_rt + 2.0e-4_rt) * dt;
    EXPECT_THROW(
        ::ext_turb::num_whole_steps(t_start, dt, "start", "dt"),
        amrex::RuntimeError);
    EXPECT_EQ(
        ::ext_turb::num_whole_steps(t_start, dt, "start", "dt", 1.0e-3_rt),
        143000);
#endif
    // Half a step off: rejected
    EXPECT_THROW(
        ::ext_turb::num_whole_steps(
            1000.5_rt * dt, dt, "start", "dt", 1.0e-3_rt),
        amrex::RuntimeError);
}

TEST_F(ExtTimeStepsTest, abort_names_inputs)
{
    std::string msg;
    try {
        ::ext_turb::num_whole_steps(
            0.00125_rt, 0.001_rt, "FastIface T0: time.fixed_dt",
            "the OpenFAST DT");
    } catch (const amrex::RuntimeError& err) {
        msg = err.what();
    }
    EXPECT_EQ(
        msg,
        "FastIface T0: time.fixed_dt = 0.00125 is not a whole multiple of "
        "the OpenFAST DT = 0.001");
}

} // namespace kynema_sgf_tests
