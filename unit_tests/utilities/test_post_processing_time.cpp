#include "ks_test_utils/MeshTest.H"
#include "src/CFDSim.H"
#include "src/utilities/PostProcessing.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

class PostProcTimeTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", -0.1_rt);
            pp.add("init_shrink", 0.1_rt);
            pp.add("cfl", 0.45_rt);
            pp.add("verbose", -1);
            pp.add("regrid_interval", -1);
            pp.add("checkpoint_interval", -1);
            pp.add("plot_interval", -1);
        }
        {
            amrex::ParmParse pp("incflo");
            pp.add("post_processing", (std::string) "fnorm");
        }
        {
            amrex::ParmParse pp("fnorm");
            pp.add("type", (std::string) "FieldNorms");
        }
    }
};

TEST_F(PostProcTimeTest, time_interval)
{
    populate_parameters();
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 0.3_rt);
        pp.add("stop_time", 5.0_rt);
        pp.add("max_step", 100);
    }
    {
        amrex::ParmParse pp("fnorm");
        pp.add("output_time_interval", 1.0_rt);
    }
    initialize_mesh();

    auto& m_sim = sim();
    kynema_sgf::PostProcessManager& post_manager = m_sim.post_manager();
    auto& time = sim().time();
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    int out_counter = 0;
    amrex::Real out_time_sum = 0.0_rt;
    int out_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.3_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        post_manager.post_advance_work();
    }

    // Read file output for time steps, times
    std::string fname = "post_processing/fnorm00000.txt";
    std::ifstream ifh(fname, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort("Cannot find file: " + fname);
    }
    int data_tstep;
    amrex::Real data_time;
    ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    while (ifh >> data_tstep) {
        ifh >> data_time;
        ++out_counter;
        out_step_sum += data_tstep;
        out_time_sum += data_time;
    }

    EXPECT_EQ(out_counter, 1 + 5);
    EXPECT_NEAR(
        out_time_sum, 0.0_rt + 1.2_rt + 2.1_rt + 3.0_rt + 4.2_rt + 5.1_rt,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
    EXPECT_EQ(out_step_sum, 4 + 7 + 10 + 14 + 17);

    // Remove file
    if (ifh.good()) {
        remove(fname.c_str());
    }
}

TEST_F(PostProcTimeTest, time_interval_perturb)
{
    populate_parameters();
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 0.3_rt);
        pp.add("stop_time", 5.0_rt);
        pp.add("max_step", 100);
    }
    {
        amrex::ParmParse pp("fnorm");
        pp.add("output_time_interval", 1.0_rt);
        // Perturb start time - like floating point error present on a restart
        pp.add("output_start_time", 1.0e-6_rt);
    }
    initialize_mesh();

    auto& m_sim = sim();
    kynema_sgf::PostProcessManager& post_manager = m_sim.post_manager();
    auto& time = sim().time();
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    int out_counter = 0;
    amrex::Real out_time_sum = 0.0_rt;
    int out_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.3_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        post_manager.post_advance_work();
    }

    // Read file output for time steps, times
    std::string fname = "post_processing/fnorm00000.txt";
    std::ifstream ifh(fname, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort("Cannot find file: " + fname);
    }
    int data_tstep;
    amrex::Real data_time;
    ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    while (ifh >> data_tstep) {
        ifh >> data_time;
        ++out_counter;
        out_step_sum += data_tstep;
        out_time_sum += data_time;
    }

    EXPECT_EQ(out_counter, 1 + 5);
    EXPECT_NEAR(
        out_time_sum, 0.0_rt + 1.2_rt + 2.1_rt + 3.0_rt + 4.2_rt + 5.1_rt,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
    EXPECT_EQ(out_step_sum, 4 + 7 + 10 + 14 + 17);

    // Remove file
    if (ifh.good()) {
        remove(fname.c_str());
    }
}

TEST_F(PostProcTimeTest, enforce_time_interval)
{
    populate_parameters();
    {
        amrex::ParmParse pp("time");
        pp.add("stop_time", 1.0_rt);
        pp.add("max_step", 20);
    }
    {
        amrex::ParmParse pp("fnorm");
        pp.add("output_time_interval", 0.45_rt);
        pp.add("enforce_output_time_dt", true);
    }
    initialize_mesh();

    auto& m_sim = sim();
    kynema_sgf::PostProcessManager& post_manager = m_sim.post_manager();
    auto& time = sim().time();
    post_manager.pre_init_actions();
    // Give a real value for dt so initial output happens
    time.delta_t() = 1.0_rt;
    post_manager.post_init_actions();

    int out_counter = 0;
    amrex::Real out_time_sum = 0.0_rt;
    int out_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.4_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        post_manager.post_advance_work();
    }

    // Read file output for time steps, times
    std::string fname = "post_processing/fnorm00000.txt";
    std::ifstream ifh(fname, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort("Cannot find file: " + fname);
    }
    int data_tstep;
    amrex::Real data_time;
    ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    while (ifh >> data_tstep) {
        ifh >> data_time;
        ++out_counter;
        out_step_sum += data_tstep;
        out_time_sum += data_time;
    }

    EXPECT_EQ(out_counter, 3);
    EXPECT_NEAR(
        out_time_sum, 0.0_rt + 0.45_rt + (2.0_rt * 0.45_rt),
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
    EXPECT_EQ(out_step_sum, 2 + 9);

    // Remove file
    if (ifh.good()) {
        remove(fname.c_str());
    }
}

TEST_F(PostProcTimeTest, output_end)
{
    populate_parameters();
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 0.3_rt);
        pp.add("stop_time", 10.0_rt);
        pp.add("max_step", 21);
    }
    {
        amrex::ParmParse pp("fnorm");
        pp.add("output_delay", 1000);
        pp.add("output_after_final_step", true);
    }
    initialize_mesh();

    auto& m_sim = sim();
    kynema_sgf::PostProcessManager& post_manager = m_sim.post_manager();
    auto& time = sim().time();
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    int out_counter = 0;
    amrex::Real out_time_sum = 0.0_rt;
    int out_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.3_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        post_manager.post_advance_work();
    }
    post_manager.final_output();

    // Read file output for time steps, times
    std::string fname = "post_processing/fnorm00000.txt";
    std::ifstream ifh(fname, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort("Cannot find file: " + fname);
    }
    int data_tstep;
    amrex::Real data_time;
    ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    while (ifh >> data_tstep) {
        ifh >> data_time;
        ++out_counter;
        out_step_sum += data_tstep;
        out_time_sum += data_time;
    }

    EXPECT_EQ(out_counter, 1);
    EXPECT_NEAR(
        out_time_sum, 6.3_rt,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
    EXPECT_EQ(out_step_sum, 21);

    // Remove file
    if (ifh.good()) {
        remove(fname.c_str());
    }
}

TEST_F(PostProcTimeTest, time_output_end)
{
    populate_parameters();
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 0.3_rt);
        pp.add("stop_time", 10.0_rt);
        pp.add("max_step", 21);
    }
    {
        amrex::ParmParse pp("fnorm");
        pp.add("output_time_interval", 1.0_rt);
        pp.add("output_time_delay", 1000.0_rt);
        pp.add("output_after_final_step", true);
    }
    initialize_mesh();

    auto& m_sim = sim();
    kynema_sgf::PostProcessManager& post_manager = m_sim.post_manager();
    auto& time = sim().time();
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    int out_counter = 0;
    amrex::Real out_time_sum = 0.0_rt;
    int out_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.3_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        post_manager.post_advance_work();
    }
    post_manager.final_output();

    // Read file output for time steps, times
    std::string fname = "post_processing/fnorm00000.txt";
    std::ifstream ifh(fname, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort("Cannot find file: " + fname);
    }
    int data_tstep;
    amrex::Real data_time;
    ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    while (ifh >> data_tstep) {
        ifh >> data_time;
        ++out_counter;
        out_step_sum += data_tstep;
        out_time_sum += data_time;
    }

    EXPECT_EQ(out_counter, 1);
    EXPECT_NEAR(
        out_time_sum, 6.3_rt,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
    EXPECT_EQ(out_step_sum, 21);

    // Remove file
    if (ifh.good()) {
        remove(fname.c_str());
    }
}

TEST_F(PostProcTimeTest, no_output_end)
{
    populate_parameters();
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 0.3_rt);
        pp.add("stop_time", 10.0_rt);
        pp.add("max_step", 21);
    }
    {
        amrex::ParmParse pp("fnorm");
        pp.add("output_delay", 1000);
    }
    initialize_mesh();

    auto& m_sim = sim();
    kynema_sgf::PostProcessManager& post_manager = m_sim.post_manager();
    auto& time = sim().time();
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    int out_counter = 0;
    amrex::Real out_time_sum = 0.0_rt;
    int out_step_sum = 0;
    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.3_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        post_manager.post_advance_work();
    }
    post_manager.final_output();

    // Read file output for time steps, times
    std::string fname = "post_processing/fnorm00000.txt";
    std::ifstream ifh(fname, std::ios::in);
    if (!ifh.good()) {
        amrex::Abort("Cannot find file: " + fname);
    }
    int data_tstep;
    amrex::Real data_time;
    ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    while (ifh >> data_tstep) {
        ifh >> data_time;
        ++out_counter;
        out_step_sum += data_tstep;
        out_time_sum += data_time;
    }

    // File should have no data and only contain the header
    EXPECT_EQ(out_counter, 0);
    EXPECT_NEAR(
        out_time_sum, 0.0_rt,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt);
    EXPECT_EQ(out_step_sum, 0);

    // Remove file
    if (ifh.good()) {
        remove(fname.c_str());
    }
}

TEST_F(PostProcTimeTest, conflict_fails)
{
    populate_parameters();
    initialize_mesh();

    auto& m_sim = sim();
    kynema_sgf::PostProcessManager& post_manager = m_sim.post_manager();
    post_manager.pre_init_actions();
    // Confirm no fail with default arguments
    post_manager.post_init_actions();
    // Confirm no fail with time interval
    {
        amrex::ParmParse pp("fnorm");
        pp.add("output_time_interval", 1.0_rt);
    }
    post_manager.post_init_actions();
    // Should fail when overspecified arguments
    {
        amrex::ParmParse pp("fnorm");
        pp.add("output_interval", 1);
    }
    EXPECT_THROW(post_manager.post_init_actions(), std::runtime_error);
    // Should fail when both invalid
    {
        amrex::ParmParse pp("fnorm");
        pp.add("output_interval", -1);
        pp.add("output_time_interval", -1.0_rt);
    }
    EXPECT_THROW(post_manager.post_init_actions(), std::runtime_error);
    // Should fail for combination of enforce and step interval
    {
        amrex::ParmParse pp("fnorm");
        pp.add("output_interval", 1);
        pp.add("enforce_output_time_dt", true);
    }
    EXPECT_THROW(post_manager.post_init_actions(), std::runtime_error);
}

} // namespace amr_wind_tests
