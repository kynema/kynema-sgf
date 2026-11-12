#include <string>
#include "ks_test_utils/AmrexTest.H"
#include "src/actuator/turbine/external/ExtTurbIface.H"
#include "AMReX_Exception.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {
namespace {

//! Time steps of an external solver in a CFD time step
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
    // In double precision 0.3 / 0.1 is 2.9999999999999996 and 0.35 / 0.05 is
    // 6.999999999999999: rounded down they give 2 and 6 steps
    EXPECT_EQ(num_steps(0.3_rt, 0.1_rt), 3);
    EXPECT_EQ(num_steps(0.0003_rt, 0.0001_rt), 3);
    EXPECT_EQ(num_steps(0.009_rt, 0.003_rt), 3);
    EXPECT_EQ(num_steps(0.35_rt, 0.05_rt), 7);
}

TEST_F(ExtTimeStepsTest, ratio_just_above_whole_number)
{
    // 0.035 / 0.005 is 7.000000000000001 in double precision
    EXPECT_EQ(num_steps(0.035_rt, 0.005_rt), 7);
    EXPECT_EQ(num_steps(0.00875_rt, 0.00125_rt), 7);
}

TEST_F(ExtTimeStepsTest, whole_ratio)
{
    EXPECT_EQ(num_steps(0.00125_rt, 0.00125_rt), 1);
    EXPECT_EQ(num_steps(0.005_rt, 0.00125_rt), 4);
    EXPECT_EQ(num_steps(0.1_rt, 0.00125_rt), 80);
    // A replay from the initial time advances the turbine by no step
    EXPECT_EQ(num_steps(0.0_rt, 0.00125_rt), 0);
}

TEST_F(ExtTimeStepsTest, rejects_ratio_that_is_not_whole)
{
    EXPECT_THROW(num_steps(0.00125_rt, 0.001_rt), amrex::RuntimeError);
    EXPECT_THROW(num_steps(0.0030001_rt, 0.001_rt), amrex::RuntimeError);
    // A solver time step larger than the CFD time step
    EXPECT_THROW(num_steps(0.0004_rt, 0.001_rt), amrex::RuntimeError);
    EXPECT_THROW(num_steps(0.0006_rt, 0.001_rt), amrex::RuntimeError);
}

TEST_F(ExtTimeStepsTest, rejects_time_step_that_is_not_positive)
{
    EXPECT_THROW(num_steps(0.005_rt, 0.0_rt), amrex::RuntimeError);
    EXPECT_THROW(num_steps(0.005_rt, -0.00125_rt), amrex::RuntimeError);
    EXPECT_THROW(num_steps(-0.005_rt, 0.00125_rt), amrex::RuntimeError);
}

TEST_F(ExtTimeStepsTest, abort_names_interval_and_time_steps)
{
    std::string msg;
    try {
        ::ext_turb::num_whole_steps(
            0.00125_rt, 0.001_rt, "FastIface: the CFD time step of T0",
            "OpenFAST time steps");
    } catch (const amrex::RuntimeError& err) {
        msg = err.what();
    }
    EXPECT_EQ(
        msg, "FastIface: the CFD time step of T0 (0.00125) is not a whole "
             "number of OpenFAST time steps (0.001)");
}

} // namespace kynema_sgf_tests
