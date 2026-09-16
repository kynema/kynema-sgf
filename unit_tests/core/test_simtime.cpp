#include "ks_test_utils/AmrexTest.H"
#include "AMReX_ParmParse.H"
#include "src/core/SimTime.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

void build_simtime_params()
{
    amrex::ParmParse pp("time");
    pp.add("stop_time", 2.0_rt);
    pp.add("max_step", 10);
    pp.add("fixed_dt", -0.1_rt);
    pp.add("init_shrink", 0.1_rt);
    pp.add("cfl", 0.45_rt);
    pp.add("verbose", -1);
    pp.add("regrid_interval", 3);
    pp.add("plot_interval", 1);
    pp.add("checkpoint_interval", 2);
}

} // namespace

//! Create unique namespace for this test fixture
class SimTimeTest : public AmrexTest
{};

TEST_F(SimTimeTest, request_stop_ends_the_run_before_max_step)
{
    build_simtime_params();
    kynema_sgf::SimTime time;
    time.parse_parameters();

    // Nothing has asked to stop, and neither max_step nor stop_time is close
    EXPECT_FALSE(time.stop_requested());
    EXPECT_TRUE(time.continue_simulation());

    time.request_stop("monitor converged");

    EXPECT_TRUE(time.stop_requested());
    EXPECT_EQ(time.stop_reason(), "monitor converged");

    // Still at the first step, well short of the max_step of 10 and the
    // stop_time of 2.0, so this is the requested stop taking effect and not
    // one of the existing end conditions
    EXPECT_EQ(time.time_index(), 0);
    EXPECT_FALSE(time.continue_simulation());
    EXPECT_FALSE(time.new_timestep());
}

namespace {

//! Advance a SimTime through a number of fixed-size steps
void advance_steps(kynema_sgf::SimTime& time, const int nsteps)
{
    for (int i = 0; i < nsteps; ++i) {
        ASSERT_TRUE(time.new_timestep());
        time.set_current_cfl(1.125_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
    }
}

} // namespace

TEST_F(SimTimeTest, requested_stop_forces_final_output_without_intervals)
{
    build_simtime_params();
    {
        // A run set up to stop on convergence and nothing else
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 0.1_rt);
        pp.add("plot_interval", -1);
        pp.add("checkpoint_interval", -1);
    }
    kynema_sgf::SimTime time;
    time.parse_parameters();
    advance_steps(time, 3);

    // The existing end-of-run rule owes nothing when no interval was set
    EXPECT_FALSE(time.write_final_plot_file());
    EXPECT_FALSE(time.write_final_checkpoint());

    // A requested stop must leave both behind regardless
    time.request_stop("monitor converged");
    EXPECT_TRUE(time.write_final_plot_file());
    EXPECT_TRUE(time.write_final_checkpoint());
}

TEST_F(SimTimeTest, requested_stop_does_not_duplicate_interval_output)
{
    // With intervals configured, the final output must be exactly what the
    // run would have written anyway: nothing on a step the interval already
    // wrote, and a last file on a step it did not. The requested stop must
    // not change that, or a stop landing on a plot step would write the
    // plot file twice
    for (int nsteps = 1; nsteps <= 6; ++nsteps) {
        build_simtime_params();
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", 0.1_rt);
            pp.add("plot_interval", 2);
            pp.add("checkpoint_interval", 3);
        }
        kynema_sgf::SimTime time;
        time.parse_parameters();
        advance_steps(time, nsteps);

        const bool plt_owed = time.write_last_plot_file();
        const bool chk_owed = time.write_last_checkpoint();
        EXPECT_EQ(plt_owed, (nsteps % 2) != 0) << "step " << nsteps;
        EXPECT_EQ(chk_owed, (nsteps % 3) != 0) << "step " << nsteps;

        time.request_stop("monitor converged");
        EXPECT_EQ(time.write_final_plot_file(), plt_owed) << "step " << nsteps;
        EXPECT_EQ(time.write_final_checkpoint(), chk_owed) << "step " << nsteps;
    }
}

TEST_F(SimTimeTest, init)
{
    build_simtime_params();
    kynema_sgf::SimTime time;
    time.parse_parameters();

    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    EXPECT_EQ(time.time_index(), 0);
    EXPECT_NEAR(time.current_time(), 0.0_rt, tol);
    EXPECT_NEAR(time.max_cfl(), 0.45_rt, tol);

    const amrex::Real cur_cfl = 0.9_rt;
    const amrex::Real dt_new =
        1.0_rt; // which comes from "2.0_rt * 0.45_rt / cur_cfl"

    // Check that the timestep size during initialization respects the shrink
    // value
    time.set_current_cfl(cur_cfl * 0.5_rt, 0.0_rt, 0.0_rt);
    EXPECT_NEAR(time.delta_t(), 0.1_rt * dt_new, tol);
    const amrex::Real first_dt = time.delta_t();

    bool stop_sim = time.new_timestep();
    EXPECT_TRUE(stop_sim);
    // Check that the timestep growth is not greater than 10% of the last
    // timestep
    time.set_current_cfl(cur_cfl * 0.5_rt, 0.0_rt, 0.0_rt);
    EXPECT_NEAR(time.delta_t(), 1.1_rt * first_dt, tol);
}

TEST_F(SimTimeTest, time_loop)
{
    build_simtime_params();
    kynema_sgf::SimTime time;
    time.parse_parameters();

    int counter = 0;
    int regrid_counter = 0;
    int plot_counter = 0;
    int chkpt_counter = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(1.125_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        ++counter;

        if (time.write_plot_file()) {
            ++plot_counter;
        }
        if (time.write_checkpoint()) {
            ++chkpt_counter;
        }
        if (time.do_regrid()) {
            ++regrid_counter;
        }
        std::cout << time.new_time() << '\n';
    }
    EXPECT_EQ(counter, 5);
    EXPECT_EQ(plot_counter, 5);
    EXPECT_EQ(chkpt_counter, 2);
    EXPECT_EQ(regrid_counter, 1);

    EXPECT_TRUE(time.write_last_checkpoint());
    EXPECT_FALSE(time.write_last_plot_file());
}

TEST_F(SimTimeTest, fixed_dt_loop)
{
    build_simtime_params();
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 0.2_rt);
    }
    kynema_sgf::SimTime time;
    time.parse_parameters();

    int counter = 0;
    int regrid_counter = 0;
    int plot_counter = 0;
    int chkpt_counter = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(2.0_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        ++counter;

        if (time.write_plot_file()) {
            ++plot_counter;
        }
        if (time.write_checkpoint()) {
            ++chkpt_counter;
        }
        if (time.do_regrid()) {
            ++regrid_counter;
        }
    }
    EXPECT_EQ(counter, 10);
    EXPECT_EQ(plot_counter, 10);
    EXPECT_EQ(chkpt_counter, 5);
    EXPECT_EQ(regrid_counter, 3);

    EXPECT_FALSE(time.write_last_checkpoint());
    EXPECT_FALSE(time.write_last_plot_file());
}

TEST_F(SimTimeTest, fixed_dt_delay)
{
    build_simtime_params();
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 0.2_rt);
        pp.add("regrid_interval", -1);
        pp.add("checkpoint_delay", 3);
        pp.add("plot_delay", 5);
    }
    kynema_sgf::SimTime time;
    time.parse_parameters();

    int plot_counter = 0;
    int chkpt_counter = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(2.0_rt, 0.0_rt, 0.0_rt);
        time.advance_time();

        if (time.write_plot_file()) {
            ++plot_counter;
        }
        if (time.write_checkpoint()) {
            ++chkpt_counter;
        }
    }
    EXPECT_EQ(plot_counter, 6);
    EXPECT_EQ(chkpt_counter, 4);
}

TEST_F(SimTimeTest, plt_chk_timeinterval_loop)
{
    build_simtime_params();
    {
        amrex::ParmParse pp("time");
        pp.add("regrid_interval", -1);
        pp.add("checkpoint_interval", -1);
        pp.add("plot_interval", -1);
        pp.add("checkpoint_time_interval", 4.0_rt);
        pp.add("plot_time_interval", 1.0_rt);
        pp.add("fixed_dt", 0.3_rt);
        pp.add("stop_time", 5.0_rt);
        pp.add("max_step", 100);
    }
    kynema_sgf::SimTime time;
    time.parse_parameters();

    int counter = 0;
    int plot_counter = 0;
    int plot_step_sum = 0;
    int chkpt_counter = 0;
    int chkpt_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.3_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        ++counter;

        if (time.write_plot_file()) {
            ++plot_counter;
            plot_step_sum += counter;
        }
        if (time.write_checkpoint()) {
            ++chkpt_counter;
            chkpt_step_sum += counter;
        }
    }
    EXPECT_EQ(plot_counter, 5);
    EXPECT_EQ(plot_step_sum, 4 + 7 + 10 + 14 + 17);
    EXPECT_EQ(chkpt_counter, 1);
    EXPECT_EQ(chkpt_step_sum, 14);
}

TEST_F(SimTimeTest, plt_chk_timeinterval_loop_perturb)
{
    build_simtime_params();
    {
        amrex::ParmParse pp("time");
        pp.add("regrid_interval", -1);
        pp.add("checkpoint_interval", -1);
        pp.add("plot_interval", -1);
        pp.add("checkpoint_time_interval", 4.0_rt);
        pp.add("plot_time_interval", 1.0_rt);
        pp.add("fixed_dt", 0.3_rt);
        pp.add("stop_time", 5.0_rt);
        pp.add("max_step", 100);
        // Perturb start time - like floating point error present on a restart
        pp.add("checkpoint_start_time", 1.0e-6_rt);
        pp.add("plot_start_time", 1.0e-6_rt);
    }
    kynema_sgf::SimTime time;
    time.parse_parameters();

    int counter = 0;
    int plot_counter = 0;
    int plot_step_sum = 0;
    int chkpt_counter = 0;
    int chkpt_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.3_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        ++counter;

        if (time.write_plot_file()) {
            ++plot_counter;
            plot_step_sum += counter;
        }
        if (time.write_checkpoint()) {
            ++chkpt_counter;
            chkpt_step_sum += counter;
        }
    }
    EXPECT_EQ(plot_counter, 5);
    EXPECT_EQ(plot_step_sum, 4 + 7 + 10 + 14 + 17);
    EXPECT_EQ(chkpt_counter, 1);
    EXPECT_EQ(chkpt_step_sum, 14);
}

TEST_F(SimTimeTest, plt_chk_timeinterval_delay)
{
    build_simtime_params();
    {
        amrex::ParmParse pp("time");
        pp.add("regrid_interval", -1);
        pp.add("checkpoint_interval", -1);
        pp.add("plot_interval", -1);
        pp.add("checkpoint_time_interval", 2.0_rt);
        pp.add("checkpoint_time_delay", 4.0_rt);
        pp.add("plot_time_interval", 1.0_rt);
        pp.add("plot_time_delay", 3.0_rt);
        pp.add("fixed_dt", 0.3_rt);
        pp.add("stop_time", 5.0_rt);
        pp.add("max_step", 100);
    }
    kynema_sgf::SimTime time;
    time.parse_parameters();

    int counter = 0;
    int plot_counter = 0;
    int plot_step_sum = 0;
    int chkpt_counter = 0;
    int chkpt_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.3_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        ++counter;

        if (time.write_plot_file()) {
            ++plot_counter;
            plot_step_sum += counter;
        }
        if (time.write_checkpoint()) {
            ++chkpt_counter;
            chkpt_step_sum += counter;
        }
    }
    EXPECT_EQ(plot_counter, 3);
    EXPECT_EQ(plot_step_sum, 10 + 14 + 17);
    EXPECT_EQ(chkpt_counter, 1);
    EXPECT_EQ(chkpt_step_sum, 14);
}

TEST_F(SimTimeTest, enforce_dt_out)
{
    // Should not change if already correct
    amrex::Real result =
        get_enforced_dt_for_output(0.1_rt, 3.9_rt, 2.0_rt, 1.0e-3_rt);
    EXPECT_NEAR(
        result, 0.1_rt, std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);

    // Should not change if short of interval
    result = get_enforced_dt_for_output(
        0.1_rt, 3.9_rt - (2.0_rt * 5.0e-4_rt), 2.0_rt, 1.0e-3_rt);
    EXPECT_NEAR(
        result, 0.1_rt, std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);

    // Should not change if starting near interval
    result = get_enforced_dt_for_output(0.1_rt, 4.0_rt, 2.0_rt, 1.0e-3_rt);
    EXPECT_NEAR(
        result, 0.1_rt, std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
    result = get_enforced_dt_for_output(
        0.1_rt, 4.0_rt - (2.0_rt * 0.99e-3_rt), 2.0_rt, 1.0e-3_rt);
    EXPECT_NEAR(
        result, 0.1_rt, std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
    // Past the tolerance, will change
    result = get_enforced_dt_for_output(
        0.1_rt, 4.0_rt - (2.0_rt * 1.01e-3_rt), 2.0_rt, 1.0e-3_rt);
    EXPECT_LT(result, 0.1_rt);

    // Shortens dt if overlapping with intervals
    result = get_enforced_dt_for_output(0.1_rt, 3.95_rt, 2.0_rt, 1.0e-3_rt);
    EXPECT_NEAR(
        result, 0.05_rt,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
}

TEST_F(SimTimeTest, enforce_timeinterval)
{
    build_simtime_params();
    {
        amrex::ParmParse pp("time");
        pp.add("regrid_interval", -1);
        pp.add("checkpoint_interval", -1);
        pp.add("plot_interval", -1);

        pp.add("plot_time_interval", 0.5_rt);
        pp.add("enforce_plot_time_dt", true);

        // Default values for tolerances
        pp.add("stop_time", 1.0_rt);
        pp.add("max_step", 10);
    }
    kynema_sgf::SimTime time;
    time.parse_parameters();

    int counter = 0;
    int plot_counter = 0;
    amrex::Real plot_time_sum = 0.0_rt;
    int plot_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.4_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        ++counter;
        if (time.write_plot_file()) {
            ++plot_counter;
            plot_time_sum += time.new_time();
            plot_step_sum += counter;
        }
    }
    EXPECT_EQ(plot_counter, 2);
    EXPECT_NEAR(
        plot_time_sum, 1.5_rt,
        std::numeric_limits<float>::epsilon() * 1.0e1_rt);
    EXPECT_EQ(plot_step_sum, 2 + 6);
}

TEST_F(SimTimeTest, enforce_timeinterval_bigtimetol)
{
    build_simtime_params();
    {
        amrex::ParmParse pp("time");
        pp.add("regrid_interval", -1);
        pp.add("checkpoint_interval", -1);
        pp.add("plot_interval", -1);

        pp.add("plot_time_interval", 0.5_rt);
        pp.add("enforce_plot_time_dt", true);

        // Choose values that could give issues
        pp.add(
            "enforce_plot_dt_reltol",
            std::numeric_limits<float>::epsilon() * 1.0e1_rt);
        pp.add("plot_time_interval_reltol", 1.0e0_rt);

        pp.add("stop_time", 1.0_rt);
        pp.add("max_step", 10);
    }
    kynema_sgf::SimTime time;
    time.parse_parameters();

    int counter = 0;
    int plot_counter = 0;
    amrex::Real plot_time_sum = 0.0_rt;
    int plot_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.4_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        ++counter;
        if (time.write_plot_file()) {
            ++plot_counter;
            plot_time_sum += time.new_time();
            plot_step_sum += counter;
        }
    }
    EXPECT_EQ(plot_counter, 2);
    EXPECT_NEAR(
        plot_time_sum, 1.5_rt,
        std::numeric_limits<float>::epsilon() * 1.0e1_rt);
    EXPECT_EQ(plot_step_sum, 2 + 6);
}

TEST_F(SimTimeTest, enforce_timeinterval_bigdttol)
{
    build_simtime_params();
    {
        amrex::ParmParse pp("time");
        pp.add("regrid_interval", -1);
        pp.add("checkpoint_interval", -1);
        pp.add("plot_interval", -1);

        pp.add("plot_time_interval", 0.5_rt);
        pp.add("enforce_plot_time_dt", true);

        // Weak enforcement of plot time interval on dt
        pp.add("enforce_plot_dt_reltol", 1.0e0_rt);
        pp.add(
            "plot_time_interval_reltol",
            std::numeric_limits<float>::epsilon() * 1.0e1_rt);

        pp.add("stop_time", 1.0_rt);
        pp.add("max_step", 10);
    }
    kynema_sgf::SimTime time;
    time.parse_parameters();

    int counter = 0;
    int plot_counter = 0;
    amrex::Real plot_time_sum = 0.0_rt;
    int plot_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.4_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        ++counter;
        if (time.write_plot_file()) {
            ++plot_counter;
            plot_time_sum += time.new_time();
            plot_step_sum += counter;
        }
    }
    // With big dt tolerance, dt never gets shortened except at the end
    // (reaching t = 1.0_rt). Ordinary plot time interval tolerance ensures that
    // plot files are still written at the first step after interval is passed.
    EXPECT_EQ(plot_counter, 2);
    EXPECT_NEAR(
        plot_time_sum, 0.8_rt + 1.0_rt,
        std::numeric_limits<float>::epsilon() * 1.0e1_rt);
    EXPECT_EQ(plot_step_sum, 2 + 3);
}

TEST_F(SimTimeTest, enforce_timeinterval_delay)
{
    build_simtime_params();
    {
        amrex::ParmParse pp("time");
        pp.add("regrid_interval", -1);
        pp.add("checkpoint_interval", -1);
        pp.add("plot_interval", -1);

        pp.add("plot_time_interval", 0.5_rt);
        pp.add("plot_time_delay", 0.9_rt);
        pp.add("enforce_plot_time_dt", true);

        // Default values for tolerances
        pp.add("stop_time", 1.0_rt);
        pp.add("max_step", 10);
    }
    kynema_sgf::SimTime time;
    time.parse_parameters();

    int counter = 0;
    int plot_counter = 0;
    amrex::Real plot_time_sum = 0.0_rt;
    int plot_step_sum = 0;
    amrex::Real time2 = 0.0_rt;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.4_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        ++counter;
        if (time.write_plot_file()) {
            ++plot_counter;
            plot_time_sum += time.new_time();
            plot_step_sum += counter;
        }
        if (counter == 2) {
            time2 = time.new_time();
        }
    }
    EXPECT_EQ(plot_counter, 1);
    EXPECT_NEAR(
        plot_time_sum, 1.0_rt,
        std::numeric_limits<float>::epsilon() * 1.0e1_rt);
    // dt should not shorten for t = 0.5_rt
    EXPECT_GT(time2, 0.5_rt);
    // leading to fewer steps
    EXPECT_EQ(plot_step_sum, 3);
}

TEST_F(SimTimeTest, enforce_chkpt_timeinterval)
{
    build_simtime_params();
    {
        amrex::ParmParse pp("time");
        pp.add("regrid_interval", -1);
        pp.add("checkpoint_interval", -1);
        pp.add("plot_interval", -1);

        pp.add("checkpoint_time_interval", 0.5_rt);
        pp.add("enforce_checkpoint_time_dt", true);

        // Default values for tolerances
        pp.add("stop_time", 1.0_rt);
        pp.add("max_step", 10);
    }
    kynema_sgf::SimTime time;
    time.parse_parameters();

    int counter = 0;
    int chkpt_counter = 0;
    amrex::Real chkpt_time_sum = 0.0_rt;
    int chkpt_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.4_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        ++counter;
        if (time.write_checkpoint()) {
            ++chkpt_counter;
            chkpt_time_sum += time.new_time();
            chkpt_step_sum += counter;
        }
    }
    EXPECT_EQ(chkpt_counter, 2);
    EXPECT_NEAR(
        chkpt_time_sum, 1.5_rt,
        std::numeric_limits<float>::epsilon() * 1.0e1_rt);
    EXPECT_EQ(chkpt_step_sum, 2 + 6);
}

} // namespace kynema_sgf_tests
