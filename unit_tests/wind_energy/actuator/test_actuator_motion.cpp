#include "src/wind_energy/actuator/motion/RigidBodyMotion.H"
#include "src/wind_energy/actuator/motion/RotorMotion.H"
#include "src/wind_energy/actuator/motion/TimeTable.H"
#include "src/utilities/constants.H"

#include "gtest/gtest.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <numbers>

#include "AMReX_ParmParse.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {
namespace {

using kynema_sgf::actuator::motion::RigidBodyMotion;
using kynema_sgf::actuator::motion::RotorMotion;
using kynema_sgf::actuator::motion::TimeTable;
using kynema_sgf::actuator::utils::ActParser;

void write_table(const std::string& filename, const std::string& contents)
{
    std::ofstream stream(filename);
    stream << contents;
}

TEST(ActuatorMotion, timetable_interpolation_derivative_and_integral)
{
    const std::string filename{"motion_scalar_table.txt"};
    write_table(filename, "Time Value\n0 0\n2 2\n");

    TimeTable table;
    table.read(filename, 1);
    EXPECT_NEAR(table.value(1.0_rt)[0], 1.0_rt, 1.0e-14_rt);
    EXPECT_NEAR(table.derivative(1.0_rt)[0], 1.0_rt, 1.0e-14_rt);
    EXPECT_NEAR(table.integral(1.0_rt)[0], 0.5_rt, 1.0e-14_rt);
    EXPECT_NEAR(table.value(3.0_rt)[0], 2.0_rt, 1.0e-14_rt);
    EXPECT_NEAR(table.integral(3.0_rt)[0], 4.0_rt, 1.0e-14_rt);

    std::remove(filename.c_str());
}

TEST(ActuatorMotion, timetable_error_extrapolation)
{
    const std::string filename{"motion_error_table.txt"};
    write_table(filename, "Time Value\n0 0\n1 1\n");
    TimeTable table;
    table.read(filename, 1, "error");
    EXPECT_THROW(static_cast<void>(table.value(2.0_rt)), amrex::RuntimeError);
    std::remove(filename.c_str());
}

TEST(ActuatorMotion, position_and_quaternion_histories)
{
    constexpr amrex::Real tol = kynema_sgf::constants::TIGHT_TOL;
    const std::string position_file{"motion_position_table.txt"};
    const std::string orientation_file{"motion_orientation_table.txt"};
    write_table(position_file, "Time X Y Z\n0 0 0 0\n1 2 0 0\n");
    write_table(orientation_file, "Time Qw Qx Qy Qz\n0 1 0 0 0\n1 0 0 0 1\n");

    amrex::ParmParse pp("MotionHistory");
    pp.add("position_timetable", position_file);
    pp.add("orientation_timetable", orientation_file);
    pp.add("orientation_format", "quaternion");
    RigidBodyMotion motion;
    motion.read_inputs(
        ActParser("UnusedMotionDefaults", "MotionHistory"),
        kynema_sgf::vs::Vector::zero(), kynema_sgf::vs::Tensor::identity());

    const auto position = motion.position(0.5_rt);
    const auto velocity = motion.translation_velocity(0.5_rt);
    const auto rotated =
        motion.orientation(0.5_rt) & kynema_sgf::vs::Vector::ihat();
    EXPECT_NEAR(position.x(), 1.0_rt, tol);
    EXPECT_NEAR(velocity.x(), 2.0_rt, tol);
    EXPECT_NEAR(rotated.x(), 0.0_rt, tol);
    EXPECT_NEAR(rotated.y(), 1.0_rt, tol);

    std::remove(position_file.c_str());
    std::remove(orientation_file.c_str());
}

TEST(ActuatorMotion, roll_pitch_yaw_and_quaternion_histories_match)
{
    constexpr amrex::Real tol = kynema_sgf::constants::TIGHT_TOL;
    const std::string rpy_file{"motion_rpy_table.txt"};
    const std::string quaternion_file{"motion_quaternion_table.txt"};
    write_table(rpy_file, "Time Roll Pitch Yaw\n0 0 0 0\n1 0 0 -180\n");
    write_table(quaternion_file, "Time Qw Qx Qy Qz\n0 1 0 0 0\n1 0 0 0 1\n");

    amrex::ParmParse rpy_pp("MotionRpy");
    rpy_pp.add("orientation_timetable", rpy_file);
    RigidBodyMotion rpy_motion;
    rpy_motion.read_inputs(
        ActParser("UnusedRpyDefaults", "MotionRpy"),
        kynema_sgf::vs::Vector::zero(), kynema_sgf::vs::Tensor::identity());

    amrex::ParmParse quaternion_pp("MotionQuaternion");
    quaternion_pp.add("orientation_timetable", quaternion_file);
    quaternion_pp.add("orientation_format", "quaternion");
    RigidBodyMotion quaternion_motion;
    quaternion_motion.read_inputs(
        ActParser("UnusedQuaternionDefaults", "MotionQuaternion"),
        kynema_sgf::vs::Vector::zero(), kynema_sgf::vs::Tensor::identity());

    for (const amrex::Real time : {0.0_rt, 0.25_rt, 0.5_rt, 1.0_rt}) {
        const auto rpy = rpy_motion.orientation(time);
        const auto quaternion = quaternion_motion.orientation(time);
        for (int n = 0; n < rpy.ncomp; ++n) {
            EXPECT_NEAR(rpy[n], quaternion[n], tol);
        }
    }

    std::remove(rpy_file.c_str());
    std::remove(quaternion_file.c_str());
}

TEST(ActuatorMotion, velocity_history_integrates_position)
{
    const std::string filename{"motion_velocity_table.txt"};
    write_table(filename, "Time Ux Uy Uz\n0 0 0 0\n2 2 0 0\n");

    amrex::ParmParse pp("MotionVelocity");
    pp.addarr("center", amrex::Vector<amrex::Real>{1.0_rt, 0.0_rt, 0.0_rt});
    pp.add("velocity_timetable", filename);
    RigidBodyMotion motion;
    motion.read_inputs(
        ActParser("UnusedVelocityDefaults", "MotionVelocity"),
        {1.0_rt, 0.0_rt, 0.0_rt}, kynema_sgf::vs::Tensor::identity());
    EXPECT_NEAR(motion.position(1.0_rt).x(), 1.5_rt, 1.0e-14_rt);
    EXPECT_NEAR(motion.translation_velocity(1.0_rt).x(), 1.0_rt, 1.0e-14_rt);

    std::remove(filename.c_str());
}

TEST(ActuatorMotion, position_and_velocity_histories_are_exclusive)
{
    const std::string position_file{"motion_exclusive_position.txt"};
    const std::string velocity_file{"motion_exclusive_velocity.txt"};
    write_table(position_file, "Time X Y Z\n0 0 0 0\n1 1 0 0\n");
    write_table(velocity_file, "Time Ux Uy Uz\n0 0 0 0\n1 1 0 0\n");

    amrex::ParmParse pp("MotionExclusive");
    pp.add("position_timetable", position_file);
    pp.add("velocity_timetable", velocity_file);
    RigidBodyMotion motion;
    EXPECT_THROW(
        motion.read_inputs(
            ActParser("UnusedExclusiveDefaults", "MotionExclusive"),
            kynema_sgf::vs::Vector::zero(), kynema_sgf::vs::Tensor::identity()),
        amrex::RuntimeError);

    std::remove(position_file.c_str());
    std::remove(velocity_file.c_str());
}

TEST(ActuatorMotion, rotor_speed_integrates_azimuth_and_holds_rate)
{
    const std::string filename{"motion_rotor_table.txt"};
    write_table(filename, "Time R1 R2\n0 0 0\n2 2 -2\n");

    amrex::ParmParse pp("RotorHistory");
    pp.add("rotor_speed_timetable", filename);
    pp.add("initial_azimuth_degrees", 90.0_rt);
    RotorMotion motion;
    motion.read_drone_inputs(
        ActParser("UnusedRotorDefaults", "RotorHistory"), 2);

    EXPECT_NEAR(motion.omega(0, 1.0_rt), 1.0_rt, 1.0e-14_rt);
    EXPECT_NEAR(
        motion.azimuth(0, 1.0_rt),
        0.5_rt + 0.5_rt * std::numbers::pi_v<amrex::Real>, 1.0e-14_rt);
    EXPECT_NEAR(
        motion.azimuth(0, 3.0_rt),
        4.0_rt + 0.5_rt * std::numbers::pi_v<amrex::Real>, 1.0e-14_rt);
    EXPECT_NEAR(
        motion.azimuth(1, 3.0_rt),
        -4.0_rt + 0.5_rt * std::numbers::pi_v<amrex::Real>, 1.0e-14_rt);

    std::remove(filename.c_str());
}

TEST(ActuatorMotion, per_rotor_initial_azimuths)
{
    amrex::ParmParse pp("RotorInitialAzimuths");
    pp.addarr("rotor_omegas", amrex::Vector<amrex::Real>{1.0_rt, -1.0_rt});
    pp.addarr(
        "initial_azimuths_degrees",
        amrex::Vector<amrex::Real>{0.0_rt, 180.0_rt});
    RotorMotion motion;
    motion.read_drone_inputs(
        ActParser("UnusedAzimuthDefaults", "RotorInitialAzimuths"), 2);
    EXPECT_NEAR(motion.azimuth(0, 0.0_rt), 0.0_rt, 1.0e-14_rt);
    EXPECT_NEAR(
        motion.azimuth(1, 0.0_rt), std::numbers::pi_v<amrex::Real>, 1.0e-14_rt);
}

} // namespace
} // namespace kynema_sgf_tests
