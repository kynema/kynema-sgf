#include "aw_test_utils/MeshTest.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "aw_test_utils/test_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

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
    amr_wind::PostProcessManager& post_manager = m_sim.post_manager();
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
        avg_val = (avg_val * (old_avg_time) + m_dt * fval) / avg_time;
        EXPECT_NEAR(max_f, avg_val, m_tol);
    }
}

TEST_F(TimeAveragingTest, phase_linear)
{
    populate_parameters();
    initialize_mesh();

    auto& m_sim = sim();
    amr_wind::PostProcessManager& post_manager = m_sim.post_manager();
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
            avg_val =
                (avg_val * (old_avg_time) + 2.0_rt * m_dt * fval) / avg_time;
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
    amr_wind::PostProcessManager& post_manager = m_sim.post_manager();
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
    amr_wind::PostProcessManager& post_manager = m_sim.post_manager();
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
    avg_val = (avg_val * (0.3_rt - 0.1_rt) + 3.0_rt * 0.1_rt) / 0.3_rt;
    avg_val = (avg_val * (0.5_rt - 0.2_rt) + 5.0_rt * 0.2_rt) / 0.5_rt;
    avg_val = (avg_val * (0.6_rt - 0.1_rt) + 6.0_rt * 0.1_rt) / 0.6_rt;
    avg_val = (avg_val * (0.8_rt - 0.2_rt) + 8.0_rt * 0.2_rt) / 0.8_rt;
    avg_val = (avg_val * (0.9_rt - 0.1_rt) + 9.0_rt * 0.1_rt) / 0.9_rt;
    EXPECT_NEAR(max_f, avg_val, m_tol);

    const amrex::Real avg_val_ideal =
        (0.2_rt * 2.0_rt + 0.1_rt * 3.0_rt + 0.2_rt * 5.0_rt + 0.1_rt * 6.0_rt +
         0.2_rt * 8.0_rt + 0.1_rt * 9.0_rt) /
        0.9_rt;

    EXPECT_NEAR(max_f, avg_val_ideal, m_tol);
}

} // namespace amr_wind_tests
