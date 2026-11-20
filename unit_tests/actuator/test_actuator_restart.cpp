#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include "ks_test_utils/AmrexTest.H"
#include "src/actuator/actuator_restart.H"
#include "src/actuator/FLLC.H"
#include "src/actuator/actuator_types.H"
#include "AMReX_REAL.H"
#include "AMReX_Utility.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {
namespace {

//! Values that need all their digits to be read back as they were
kynema_sgf::actuator::VecList awkward_values()
{
    namespace vs = kynema_sgf::vs;
    kynema_sgf::actuator::VecList values(3);
    values[0] = vs::Vector{1.0_rt / 3.0_rt, -2.0_rt / 7.0_rt, 1.0e-17_rt};
    values[1] = vs::Vector{-314159.26535897932_rt, 0.1_rt, 0.0_rt};
    values[2] = vs::Vector{6.02214076e23_rt, -1.0e-300_rt, 4.19_rt};
    return values;
}

/** FLLC data as fllc_init leaves them: ``npts`` force, ``npts_vel``
 *  velocity and ``npts_n`` non-uniform points; values depend on ``seed``
 */
kynema_sgf::actuator::FLLCData make_fllc_data(
    const int npts,
    const int npts_vel,
    const int npts_n,
    const amrex::Real seed)
{
    namespace vs = kynema_sgf::vs;
    kynema_sgf::actuator::FLLCData data;
    data.eps_dr = 3.0_rt;
    data.nonuniform = true;
    data.mismatched_points = (npts != npts_vel);
    data.correction_velocity.resize(npts);
    data.optimal_epsilon.resize(npts);
    data.r.resize(npts);
    data.dr.resize(npts);
    data.span_distance_force.resize(npts);
    data.span_distance_vel.resize(npts_vel);
    for (int i = 0; i < npts; ++i) {
        const auto x = static_cast<amrex::Real>(i + 1);
        data.correction_velocity[i] =
            vs::Vector{seed / x, -seed * x / 7.0_rt, seed * 1.0e-9_rt};
        // The optimal epsilon comes from the chord, not from the positions
        data.optimal_epsilon[i] = 0.25_rt / x;
        data.r[i] = seed * x / 3.0_rt;
        data.dr[i] = seed / 3.0_rt;
        data.span_distance_force[i] = seed * x / 3.0_rt;
    }
    for (int i = 0; i < npts_vel; ++i) {
        data.span_distance_vel[i] = seed * static_cast<amrex::Real>(i) / 7.0_rt;
    }
    data.nonuniform_r.resize(npts_n);
    data.nonuniform_dr.resize(npts_n);
    data.nonuniform_vel_rel.resize(npts_n);
    data.nonuniform_optimal_epsilon.resize(npts_n);
    data.nonuniform_lift.resize(npts_n);
    for (int i = 0; i < npts_n; ++i) {
        data.nonuniform_r[i] = seed * static_cast<amrex::Real>(i) / 11.0_rt;
        data.nonuniform_dr[i] = seed / 11.0_rt;
    }
    data.initialized = true;
    return data;
}

} // namespace

TEST(ActuatorRestart, fllc_state_is_read_back_as_written)
{
    // The state written at the checkpoint: 5 force points, 6 velocity
    // points, 9 non-uniform points
    const auto written = make_fllc_data(5, 6, 9, 1.0_rt);
    std::stringstream ss;
    ss << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10);
    kynema_sgf::actuator::fllc_write_state(ss, written);

    // Restarted-run data: distances differ in the last digits, one fewer
    // non-uniform point
    auto data = make_fllc_data(5, 6, 8, 1.0_rt + 1.0e-7_rt);
    data.mismatched_points = false;
    for (auto& du : data.correction_velocity) {
        du = kynema_sgf::vs::Vector::zero();
    }
    ASSERT_TRUE(kynema_sgf::actuator::fllc_read_state(ss, data));

    EXPECT_TRUE(data.mismatched_points);
    ASSERT_EQ(data.nonuniform_r.size(), written.nonuniform_r.size());
    EXPECT_EQ(data.nonuniform_dr.size(), written.nonuniform_r.size());
    EXPECT_EQ(data.nonuniform_vel_rel.size(), written.nonuniform_r.size());
    EXPECT_EQ(
        data.nonuniform_optimal_epsilon.size(), written.nonuniform_r.size());
    EXPECT_EQ(data.nonuniform_lift.size(), written.nonuniform_r.size());
    for (int i = 0; i < static_cast<int>(written.r.size()); ++i) {
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            EXPECT_EQ(
                data.correction_velocity[i][n],
                written.correction_velocity[i][n]);
        }
        EXPECT_EQ(data.r[i], written.r[i]);
        EXPECT_EQ(data.dr[i], written.dr[i]);
        EXPECT_EQ(data.span_distance_force[i], written.span_distance_force[i]);
    }
    for (int i = 0; i < static_cast<int>(written.span_distance_vel.size());
         ++i) {
        EXPECT_EQ(data.span_distance_vel[i], written.span_distance_vel[i]);
    }
    for (int i = 0; i < static_cast<int>(written.nonuniform_r.size()); ++i) {
        EXPECT_EQ(data.nonuniform_r[i], written.nonuniform_r[i]);
        EXPECT_EQ(data.nonuniform_dr[i], written.nonuniform_dr[i]);
    }
}

TEST(ActuatorRestart, fllc_state_that_does_not_fit_is_refused)
{
    const auto written = make_fllc_data(5, 6, 9, 1.0_rt);
    std::stringstream ss;
    ss << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10);
    kynema_sgf::actuator::fllc_write_state(ss, written);
    const std::string state = ss.str();

    // Restarted-run distances 1 % off (within turning/bending); each case
    // below is refused for the reason given

    // Another number of force points
    {
        auto data = make_fllc_data(6, 6, 9, 1.01_rt);
        std::stringstream is(state);
        EXPECT_FALSE(kynema_sgf::actuator::fllc_read_state(is, data));
        EXPECT_EQ(data.r[0], 1.01_rt / 3.0_rt);
    }
    // Another number of velocity points
    {
        auto data = make_fllc_data(5, 5, 9, 1.01_rt);
        std::stringstream is(state);
        EXPECT_FALSE(kynema_sgf::actuator::fllc_read_state(is, data));
        EXPECT_EQ(data.r[0], 1.01_rt / 3.0_rt);
    }
    // Another epsilon_chord or chord
    {
        auto data = make_fllc_data(5, 6, 9, 1.01_rt);
        data.optimal_epsilon[2] *= 1.01_rt;
        std::stringstream is(state);
        EXPECT_FALSE(kynema_sgf::actuator::fllc_read_state(is, data));
        EXPECT_EQ(data.r[0], 1.01_rt / 3.0_rt);
    }
    // Another fllc_epsilon_dr_ratio
    {
        auto data = make_fllc_data(5, 6, 9, 1.01_rt);
        data.eps_dr = 2.0_rt;
        std::stringstream is(state);
        EXPECT_FALSE(kynema_sgf::actuator::fllc_read_state(is, data));
        EXPECT_EQ(data.r[0], 1.01_rt / 3.0_rt);
    }
    // The same line and the same inputs
    {
        auto data = make_fllc_data(5, 6, 9, 1.01_rt);
        std::stringstream is(state);
        EXPECT_TRUE(kynema_sgf::actuator::fllc_read_state(is, data));
        EXPECT_EQ(data.r[0], 1.0_rt / 3.0_rt);
    }
}

TEST(ActuatorRestart, fllc_state_of_other_positions_is_refused)
{
    const auto written = make_fllc_data(5, 6, 9, 1.0_rt);
    std::stringstream ss;
    ss << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10);
    kynema_sgf::actuator::fllc_write_state(ss, written);
    const std::string state = ss.str();

    // The same numbers of points and the same chords on a line that is 5 %
    // longer
    {
        auto data = make_fllc_data(5, 6, 9, 1.05_rt);
        std::stringstream is(state);
        EXPECT_FALSE(kynema_sgf::actuator::fllc_read_state(is, data));
        EXPECT_EQ(data.r[0], 1.05_rt / 3.0_rt);
    }
    // The same force points, and a velocity point at another place
    {
        auto data = make_fllc_data(5, 6, 9, 1.0_rt);
        data.span_distance_vel[3] += 0.1_rt;
        std::stringstream is(state);
        EXPECT_FALSE(kynema_sgf::actuator::fllc_read_state(is, data));
        EXPECT_EQ(data.span_distance_vel[3], (3.0_rt / 7.0_rt) + 0.1_rt);
    }
    // Distances 0.5 % away from those of the state: the same line
    {
        auto data = make_fllc_data(5, 6, 9, 1.005_rt);
        std::stringstream is(state);
        EXPECT_TRUE(kynema_sgf::actuator::fllc_read_state(is, data));
        EXPECT_EQ(data.r[0], 1.0_rt / 3.0_rt);
    }
}

TEST(ActuatorRestart, list_is_read_back_as_written)
{
    namespace restart = kynema_sgf::actuator::restart;
    const auto values = awkward_values();
    std::stringstream ss;
    ss << std::setprecision(std::numeric_limits<amrex::Real>::max_digits10);
    restart::write_list(ss, "correction_velocity", values);

    kynema_sgf::actuator::VecList back(values.size());
    ASSERT_TRUE(restart::read_list(ss, "correction_velocity", back));
    for (int i = 0; i < static_cast<int>(values.size()); ++i) {
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            EXPECT_EQ(back[i][n], values[i][n]) << "entry " << i << ", " << n;
        }
    }
}

TEST(ActuatorRestart, list_of_another_name_or_size_is_refused)
{
    namespace vs = kynema_sgf::vs;
    namespace restart = kynema_sgf::actuator::restart;
    const auto values = awkward_values();
    std::stringstream ss;
    restart::write_list(ss, "correction_velocity", values);

    // A state written for another number of points
    kynema_sgf::actuator::VecList longer(
        values.size() + 1, vs::Vector{7.0_rt, 7.0_rt, 7.0_rt});
    ss.seekg(0);
    EXPECT_FALSE(restart::read_list(ss, "correction_velocity", longer));
    EXPECT_EQ(longer[0][0], 7.0_rt);

    // Another list than the one asked for
    kynema_sgf::actuator::VecList same_size(
        values.size(), vs::Vector{7.0_rt, 7.0_rt, 7.0_rt});
    ss.clear();
    ss.seekg(0);
    EXPECT_FALSE(restart::read_list(ss, "fluid_velocity", same_size));
    EXPECT_EQ(same_size[0][0], 7.0_rt);

    // A list that ends early
    std::stringstream cut("correction_velocity 3\n1.0 2.0 3.0\n4.0 5.0\n");
    EXPECT_FALSE(restart::read_list(cut, "correction_velocity", same_size));
    EXPECT_EQ(same_size[0][0], 7.0_rt);
}

TEST(ActuatorRestart, file_is_for_one_actuator_and_one_time)
{
    namespace restart = kynema_sgf::actuator::restart;
    const std::string chkname = "test_actuator_restart_chk";
    amrex::UtilCreateDirectory(chkname, 0755);
    const amrex::Real time = 1.0_rt / 3.0_rt;
    const auto values = awkward_values();
    {
        auto os = restart::create_file(chkname, "T0", time);
        restart::write_value(os, "num_blades", 3);
        restart::write_list(os, "correction_velocity", values);
    }

    std::ifstream is;
    ASSERT_TRUE(restart::open_file(is, chkname, "T0", time));
    int num_blades = 0;
    EXPECT_FALSE(restart::read_value(is, "num_points", num_blades));
    is.clear();
    is.close();

    std::ifstream is_again;
    ASSERT_TRUE(restart::open_file(is_again, chkname, "T0", time));
    ASSERT_TRUE(restart::read_value(is_again, "num_blades", num_blades));
    EXPECT_EQ(num_blades, 3);
    kynema_sgf::actuator::VecList back(values.size());
    ASSERT_TRUE(restart::read_list(is_again, "correction_velocity", back));
    for (int i = 0; i < static_cast<int>(values.size()); ++i) {
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            EXPECT_EQ(back[i][n], values[i][n]) << "entry " << i << ", " << n;
        }
    }

    // The state of another time, of another actuator, of no checkpoint
    std::ifstream is_time;
    EXPECT_FALSE(restart::open_file(is_time, chkname, "T0", 2.0_rt * time));
    std::ifstream is_label;
    EXPECT_FALSE(restart::open_file(is_label, chkname, "T1", time));
    std::ifstream is_none;
    EXPECT_FALSE(restart::open_file(is_none, "no_such_checkpoint", "T0", time));
}

TEST(ActuatorRestart, checkpoint_has_file_only_of_actuator_that_wrote_one)
{
    namespace restart = kynema_sgf::actuator::restart;
    const std::string chkname = "test_actuator_restart_has_file_chk";
    amrex::UtilCreateDirectory(chkname, 0755);
    const amrex::Real time = 0.1_rt;
    {
        auto os = restart::create_file(chkname, "T0", time);
        restart::write_value(os, "num_blades", 3);
    }

    // At the first restart force call: T0 has a state; T1 none (FLLC not
    // started); nor do other times or checkpoints without actuator data.
    EXPECT_TRUE(restart::has_file(chkname, "T0", time));
    EXPECT_FALSE(restart::has_file(chkname, "T1", time));
    EXPECT_FALSE(restart::has_file(chkname, "T0", time + 0.00125_rt));
    EXPECT_FALSE(restart::has_file("no_such_checkpoint", "T0", time));
}

} // namespace kynema_sgf_tests
