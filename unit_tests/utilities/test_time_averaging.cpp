#include <algorithm>
#include <cmath>

#include "ks_test_utils/MeshTest.H"
#include "src/CFDSim.H"
#include "src/utilities/PostProcessing.H"
#include "ks_test_utils/test_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

class TimeAveragingTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", m_dt);
        }
        {
            amrex::ParmParse pp("incflo");
            pp.add("post_processing", (std::string) "tavg");
        }
        {
            amrex::ParmParse pp("tavg");
            pp.add("type", (std::string) "TimeAveraging");
            pp.add("labels", (std::string) "means");
            pp.add("averaging_window", (amrex::Real)m_fwidth);
            pp.add("averaging_time_interval", 2.0_rt * m_dt);
        }
        {
            amrex::ParmParse pp("tavg.means");
            pp.add("fields", (std::string) "temperature");
            pp.add("averaging_type", (std::string) "ReAveraging");
        }
    }

    //! Average every step with a given window and number of time steps
    void use_window(const amrex::Real window, const int max_step)
    {
        populate_parameters();
        {
            amrex::ParmParse pp("time");
            pp.add("stop_time", 100.0_rt);
            pp.add("max_step", max_step);
        }
        {
            amrex::ParmParse pp("tavg");
            pp.add("averaging_window", window);
            pp.add("averaging_time_interval", -1.0_rt);
        }
    }

    const std::string m_name = "temperature_mean_tavg";
    const amrex::Real m_fwidth = 2.0_rt;
    const amrex::Real m_dt = 0.1_rt;
    const amrex::Real m_tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e8_rt;
};

TEST_F(TimeAveragingTest, every_step)
{
    populate_parameters();
    {
        amrex::ParmParse pp("tavg");
        pp.add("averaging_time_interval", -1.0_rt);
    }
    initialize_mesh();

    auto& m_sim = sim();
    kynema_sgf::PostProcessManager& post_manager = m_sim.post_manager();
    auto& time = sim().time();
    auto& temp = sim().repo().declare_cc_field("temperature", 1, 1, 1);
    temp.setVal(0.0_rt);
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    const auto& f_avg = sim().repo().get_field(m_name);

    amrex::Real avg_val{0.0_rt};
    while (time.new_timestep()) {
        // Give field a linear profile with time
        const amrex::Real fval = 10.0_rt * time.current_time();
        temp.setVal(fval);
        time.advance_time();
        post_manager.post_advance_work();
        // Get field values
        const auto max_f = utils::field_max(f_avg);
        const auto min_f = utils::field_min(f_avg);
        // Confirm mean fields are uniform
        EXPECT_NEAR(max_f, min_f, m_tol);
        // Confirm value is as expected
        const amrex::Real avg_time = amrex::max<amrex::Real>(
            amrex::min<amrex::Real>(time.new_time(), m_fwidth), m_dt);
        const amrex::Real old_avg_time =
            amrex::max<amrex::Real>(avg_time - m_dt, 0.0_rt);
        avg_val = ((avg_val * (old_avg_time)) + (m_dt * fval)) / avg_time;
        EXPECT_NEAR(max_f, avg_val, m_tol);
    }
}

TEST_F(TimeAveragingTest, phase_linear)
{
    populate_parameters();
    initialize_mesh();

    auto& m_sim = sim();
    kynema_sgf::PostProcessManager& post_manager = m_sim.post_manager();
    auto& time = sim().time();
    auto& temp = sim().repo().declare_cc_field("temperature", 1, 1, 1);
    temp.setVal(0.0);
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    const auto& f_avg = sim().repo().get_field(m_name);

    amrex::Real avg_val{0.0};
    int step_count{0};
    while (time.new_timestep()) {
        ++step_count;
        // Give field a linear profile with time
        const amrex::Real fval = 10.0_rt * time.current_time();
        temp.setVal(fval);
        time.advance_time();
        post_manager.post_advance_work();
        // Get field values
        const auto max_f = utils::field_max(f_avg);
        const auto min_f = utils::field_min(f_avg);
        // Confirm mean fields are uniform
        EXPECT_NEAR(max_f, min_f, m_tol);
        // Confirm value is as expected
        if (step_count % 2 == 0) {
            const amrex::Real avg_time = amrex::max<amrex::Real>(
                amrex::min<amrex::Real>(time.new_time(), m_fwidth),
                2.0_rt * m_dt);
            const amrex::Real old_avg_time =
                amrex::max<amrex::Real>(avg_time - (2.0_rt * m_dt), 0.0_rt);
            avg_val = ((avg_val * (old_avg_time)) + (2.0_rt * m_dt * fval)) /
                      avg_time;
            EXPECT_NEAR(max_f, avg_val, m_tol);
        }
    }
}

TEST_F(TimeAveragingTest, phase_stairstep_offset)
{
    populate_parameters();
    {
        amrex::ParmParse pp("tavg");
        pp.add("averaging_start_time", 0.1_rt);
    }
    initialize_mesh();

    auto& m_sim = sim();
    kynema_sgf::PostProcessManager& post_manager = m_sim.post_manager();
    auto& time = sim().time();
    auto& temp = sim().repo().declare_cc_field("temperature", 1, 1, 1);
    temp.setVal(0.0_rt);
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    const auto& f_avg = sim().repo().get_field(m_name);

    int step_count{0};
    while (time.new_timestep()) {
        ++step_count;
        // Give field a linear profile with time
        const amrex::Real fval = (step_count % 2) * 5.0_rt;
        temp.setVal(fval);
        time.advance_time();
        post_manager.post_advance_work();
        // Get field values
        const auto max_f = utils::field_max(f_avg);
        const auto min_f = utils::field_min(f_avg);
        // Confirm mean fields are uniform
        EXPECT_NEAR(max_f, min_f, m_tol);
        // Confirm value is as expected
        if (step_count % 2 != 0) {
            EXPECT_NEAR(max_f, 5.0_rt, m_tol);
        }
    }
}

TEST_F(TimeAveragingTest, mismatch_time_interval)
{
    populate_parameters();
    {
        amrex::ParmParse pp("tavg");
        pp.add("averaging_time_interval", 0.15_rt);
    }
    initialize_mesh();

    auto& m_sim = sim();
    kynema_sgf::PostProcessManager& post_manager = m_sim.post_manager();
    auto& time = sim().time();
    auto& temp = sim().repo().declare_cc_field("temperature", 1, 1, 1);
    temp.setVal(0.0);
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    const auto& f_avg = sim().repo().get_field(m_name);

    while (time.new_timestep()) {
        time.advance_time();
        // Give field a linear profile with time
        const amrex::Real fval = 10.0_rt * time.new_time();
        temp.setVal(fval);
        post_manager.post_advance_work();
    }

    // Get field values
    const auto max_f = utils::field_max(f_avg);
    const auto min_f = utils::field_min(f_avg);
    // Confirm mean fields are uniform
    EXPECT_NEAR(max_f, min_f, m_tol);
    // Times when averaging should occur:
    // 0.2_rt, 0.3_rt, 0.5_rt, 0.6_rt, 0.8_rt, 0.9_rt
    amrex::Real avg_val = 0.2_rt * 2.0_rt / 0.2_rt;
    avg_val = ((avg_val * (0.3_rt - 0.1_rt)) + (3.0_rt * 0.1_rt)) / 0.3_rt;
    avg_val = ((avg_val * (0.5_rt - 0.2_rt)) + (5.0_rt * 0.2_rt)) / 0.5_rt;
    avg_val = ((avg_val * (0.6_rt - 0.1_rt)) + (6.0_rt * 0.1_rt)) / 0.6_rt;
    avg_val = ((avg_val * (0.8_rt - 0.2_rt)) + (8.0_rt * 0.2_rt)) / 0.8_rt;
    avg_val = ((avg_val * (0.9_rt - 0.1_rt)) + (9.0_rt * 0.1_rt)) / 0.9_rt;
    EXPECT_NEAR(max_f, avg_val, m_tol);

    const amrex::Real avg_val_ideal =
        ((0.2_rt * 2.0_rt) + (0.1_rt * 3.0_rt) + (0.2_rt * 5.0_rt) +
         (0.1_rt * 6.0_rt) + (0.2_rt * 8.0_rt) + (0.1_rt * 9.0_rt)) /
        0.9_rt;

    EXPECT_NEAR(max_f, avg_val_ideal, m_tol);
}

//! Before the window is filled the average is the arithmetic time mean
TEST_F(TimeAveragingTest, window_not_filled_is_arithmetic_mean)
{
    const amrex::Real window = 1.0_rt;
    const int nsteps = 10; // nsteps * dt == window
    use_window(window, nsteps);
    initialize_mesh();

    auto& m_sim = sim();
    kynema_sgf::PostProcessManager& post_manager = m_sim.post_manager();
    auto& time = sim().time();
    auto& temp = sim().repo().declare_cc_field("temperature", 1, 1, 1);
    temp.setVal(0.0_rt);
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    const auto& f_avg = sim().repo().get_field(m_name);

    int step_count{0};
    while (time.new_timestep()) {
        ++step_count;
        // Sample n at step n, so the arithmetic mean of 1..n is (n + 1) / 2
        temp.setVal(static_cast<amrex::Real>(step_count));
        time.advance_time();
        post_manager.post_advance_work();

        const auto max_f = utils::field_max(f_avg);
        const auto min_f = utils::field_min(f_avg);
        EXPECT_NEAR(max_f, min_f, m_tol);
        EXPECT_NEAR(max_f, 0.5_rt * (step_count + 1), m_tol);
    }
    EXPECT_EQ(step_count, nsteps);
}

//! A constant field must average to that constant before and after the window
//! is exceeded (i.e., no 1 - exp(-t/W) ramp-up)
TEST_F(TimeAveragingTest, window_constant_field)
{
    const amrex::Real window = 1.0_rt;
    const int nsteps = 30; // 3 windows
    use_window(window, nsteps);
    initialize_mesh();

    auto& m_sim = sim();
    kynema_sgf::PostProcessManager& post_manager = m_sim.post_manager();
    auto& time = sim().time();
    auto& temp = sim().repo().declare_cc_field("temperature", 1, 1, 1);
    temp.setVal(1.0_rt);
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    const auto& f_avg = sim().repo().get_field(m_name);

    int step_count{0};
    while (time.new_timestep()) {
        ++step_count;
        temp.setVal(1.0_rt);
        time.advance_time();
        post_manager.post_advance_work();

        EXPECT_NEAR(utils::field_max(f_avg), 1.0_rt, m_tol);
        EXPECT_NEAR(utils::field_min(f_avg), 1.0_rt, m_tol);
    }
    EXPECT_EQ(step_count, nsteps);
}

//! Once the window is exceeded the average is an exponentially weighted
//! moving average, not an arithmetic mean over the last window:
//!   avg_new = avg_old * (1 - dt / W) + f * dt / W
//! For a step from 0 to 1, after k steps avg = 1 - (1 - dt / W)^k. After one
//! window (k = W / dt) the average is ~1 - 1/e, not 1.
TEST_F(TimeAveragingTest, window_step_response_is_exponential)
{
    const amrex::Real window = 1.0_rt;
    const int nsteps = 50;
    const int step_up = 20; // elapsed time is already > window at the step
    const int nwindow = static_cast<int>(std::round(window / m_dt));
    use_window(window, nsteps);
    initialize_mesh();

    auto& m_sim = sim();
    kynema_sgf::PostProcessManager& post_manager = m_sim.post_manager();
    auto& time = sim().time();
    auto& temp = sim().repo().declare_cc_field("temperature", 1, 1, 1);
    temp.setVal(0.0_rt);
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    const auto& f_avg = sim().repo().get_field(m_name);

    const amrex::Real decay = 1.0_rt - (m_dt / window);
    int step_count{0};
    while (time.new_timestep()) {
        ++step_count;
        temp.setVal(step_count > step_up ? 1.0_rt : 0.0_rt);
        time.advance_time();
        post_manager.post_advance_work();

        const auto max_f = utils::field_max(f_avg);
        const auto min_f = utils::field_min(f_avg);
        EXPECT_NEAR(max_f, min_f, m_tol);

        const int k = std::max(step_count - step_up, 0);
        EXPECT_NEAR(max_f, 1.0_rt - std::pow(decay, k), m_tol);

        if (k == nwindow) {
            // An arithmetic mean over the last window would be exactly 1.
            // Here ~35% of the weight still comes from before the step.
            EXPECT_LT(max_f, 0.7_rt);
            EXPECT_GT(max_f, 0.6_rt);
        }
    }
    EXPECT_EQ(step_count, nsteps);
}

} // namespace kynema_sgf_tests
