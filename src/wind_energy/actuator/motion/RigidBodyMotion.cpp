#include "src/wind_energy/actuator/motion/RigidBodyMotion.H"

#include "src/core/vs/quaternion.H"
#include "src/wind_energy/actuator/sector/actuator_sector_ops.H"

#include <algorithm>
#include <cmath>
#include <limits>

#include "AMReX.H"

using namespace amrex::literals;

namespace kynema_sgf::actuator::motion {
namespace {

vs::Vector vector3(const RealList& values)
{
    return {values[0], values[1], values[2]};
}

vs::Tensor integrate_global_angular_velocity(
    const TimeTable& table,
    const vs::Tensor& initial_orientation,
    const amrex::Real time)
{
    const auto& times = table.times();
    vs::Tensor result = initial_orientation;
    amrex::Real start = times.front();
    const amrex::Real direction = (time >= start) ? 1.0_rt : -1.0_rt;
    // Compose global-frame rotations one table interval at a time. Midpoint
    // sampling gives a second-order update for a varying angular velocity.
    while (direction * (time - start) > 0.0_rt) {
        amrex::Real end = time;
        if (direction > 0.0_rt) {
            const auto next =
                std::upper_bound(times.begin(), times.end(), start);
            if (next != times.end()) {
                end = std::min(time, *next);
            }
        } else {
            end = time;
        }
        const amrex::Real midpoint = 0.5_rt * (start + end);
        const auto omega = vector3(table.value(midpoint));
        result =
            sector::rotation_matrix_from_vector(omega * (end - start)) & result;
        start = end;
    }
    return result;
}

} // namespace

void RigidBodyMotion::read_inputs(
    const utils::ActParser& pp,
    const vs::Vector& initial_position,
    const vs::Tensor& initial_orientation,
    const std::string& angular_velocity_key,
    const vs::Vector& default_angular_velocity)
{
    m_initial_position = initial_position;
    m_initial_orientation = initial_orientation;
    m_constant_angular_velocity = default_angular_velocity;
    std::string extrapolation{"hold"};
    pp.query("timetable_extrapolation", extrapolation);

    // Translation may be prescribed by position or velocity, but never both.
    const bool has_position = pp.contains("position_timetable");
    const bool has_velocity_table = pp.contains("velocity_timetable");
    const bool has_velocity = pp.contains("translation_velocity");
    if (static_cast<int>(has_position) + static_cast<int>(has_velocity_table) +
            static_cast<int>(has_velocity) >
        1) {
        amrex::Abort(
            "Specify only one of position_timetable, velocity_timetable, and "
            "translation_velocity");
    }
    if (has_position && pp.contains("center")) {
        amrex::Abort("center cannot be combined with position_timetable");
    }
    if (has_velocity_table && !pp.contains("center")) {
        amrex::Abort("center is required with velocity_timetable");
    }
    if (has_position) {
        std::string filename;
        pp.get("position_timetable", filename);
        m_position_table.read(filename, 3, extrapolation);
    } else if (has_velocity_table) {
        std::string filename;
        pp.get("velocity_timetable", filename);
        m_velocity_table.read(filename, 3, extrapolation);
    } else {
        pp.query("translation_velocity", m_constant_velocity);
    }

    // Rotation follows the same exclusive-source rule as translation.
    const bool has_orientation = pp.contains("orientation_timetable");
    const bool has_angular_table = pp.contains("angular_velocity_timetable");
    const bool has_angular = pp.contains(angular_velocity_key);
    if (static_cast<int>(has_orientation) +
            static_cast<int>(has_angular_table) +
            static_cast<int>(has_angular) >
        1) {
        amrex::Abort(
            "Specify only one of orientation_timetable, "
            "angular_velocity_timetable, and constant angular velocity");
    }
    std::string angular_frame{"global"};
    pp.query("angular_velocity_frame", angular_frame);
    if (amrex::toLower(angular_frame) != "global") {
        amrex::Abort(
            "The first motion implementation supports only global "
            "angular_velocity_frame");
    }
    if (has_orientation) {
        std::string filename;
        pp.get("orientation_timetable", filename);
        std::string format{"roll_pitch_yaw"};
        pp.query("orientation_format", format);
        format = amrex::toLower(format);
        if (format == "roll_pitch_yaw") {
            m_orientation_format = OrientationFormat::RollPitchYaw;
            m_orientation_table.read(filename, 3, extrapolation);
        } else if (format == "quaternion") {
            m_orientation_format = OrientationFormat::Quaternion;
            m_orientation_table.read(filename, 4, extrapolation);
        } else {
            amrex::Abort(
                "orientation_format must be roll_pitch_yaw or quaternion");
        }
    } else if (has_angular_table) {
        std::string filename;
        pp.get("angular_velocity_timetable", filename);
        m_angular_velocity_table.read(filename, 3, extrapolation);
    } else {
        pp.query(angular_velocity_key, m_constant_angular_velocity);
    }
}

vs::Vector RigidBodyMotion::position(const amrex::Real time) const
{
    if (!m_position_table.empty()) {
        return vector3(m_position_table.value(time));
    }
    if (!m_velocity_table.empty()) {
        return m_initial_position + vector3(m_velocity_table.integral(time));
    }
    return m_initial_position + m_constant_velocity * time;
}

vs::Vector RigidBodyMotion::translation_velocity(const amrex::Real time) const
{
    if (!m_position_table.empty()) {
        return vector3(m_position_table.derivative(time));
    }
    if (!m_velocity_table.empty()) {
        return vector3(m_velocity_table.value(time));
    }
    return m_constant_velocity;
}

vs::Quaternion RigidBodyMotion::orientation_quaternion(const int row) const
{
    const auto values = m_orientation_table.row(row);
    if (m_orientation_format == OrientationFormat::Quaternion) {
        return vs::normalized({values[0], values[1], values[2], values[3]});
    }
    return vs::from_roll_pitch_yaw({values[0], values[1], values[2]});
}

vs::Tensor RigidBodyMotion::orientation(const amrex::Real time) const
{
    if (!m_orientation_table.empty()) {
        const auto& times = m_orientation_table.times();
        if (time <= times.front()) {
            static_cast<void>(m_orientation_table.value(time));
            return vs::tensor(orientation_quaternion(0));
        }
        if (time >= times.back()) {
            static_cast<void>(m_orientation_table.value(time));
            return vs::tensor(
                orientation_quaternion(static_cast<int>(times.size()) - 1));
        }
        const int upper = static_cast<int>(
            std::upper_bound(times.begin(), times.end(), time) - times.begin());
        const int lower = upper - 1;
        const amrex::Real fraction =
            (time - times[lower]) / (times[upper] - times[lower]);
        return vs::tensor(
            vs::slerp(
                orientation_quaternion(lower), orientation_quaternion(upper),
                fraction));
    }
    if (!m_angular_velocity_table.empty()) {
        return integrate_global_angular_velocity(
            m_angular_velocity_table, m_initial_orientation, time);
    }
    return sector::rotation_matrix_from_vector(
               m_constant_angular_velocity * time) &
           m_initial_orientation;
}

vs::Vector RigidBodyMotion::angular_velocity(const amrex::Real time) const
{
    if (!m_orientation_table.empty()) {
        const auto& times = m_orientation_table.times();
        if (times.size() == 1 || time < times.front() || time > times.back()) {
            static_cast<void>(m_orientation_table.value(time));
            return vs::Vector::zero();
        }
        int upper = static_cast<int>(
            std::upper_bound(times.begin(), times.end(), time) - times.begin());
        upper = std::min(upper, static_cast<int>(times.size()) - 1);
        const int lower = upper - 1;
        auto a = orientation_quaternion(lower);
        auto b = orientation_quaternion(upper);
        if (vs::dot(a, b) < 0.0_rt) {
            b = {-b.w, -b.x, -b.y, -b.z};
        }
        const auto relative = vs::normalized(b * vs::conjugate(a));
        // Convert the relative quaternion over this interval into a constant
        // body-frame angular velocity, then express it in the CFD frame.
        const amrex::Real half_angle =
            std::acos(std::clamp(relative.w, -1.0_rt, 1.0_rt));
        const amrex::Real sine = std::sin(half_angle);
        if (std::abs(sine) <= std::numeric_limits<amrex::Real>::epsilon()) {
            return vs::Vector::zero();
        }
        const amrex::Real scale =
            2.0_rt * half_angle / (sine * (times[upper] - times[lower]));
        const vs::Vector body_omega{
            relative.x * scale, relative.y * scale, relative.z * scale};
        return orientation(time) & body_omega;
    }
    if (!m_angular_velocity_table.empty()) {
        return vector3(m_angular_velocity_table.value(time));
    }
    return m_constant_angular_velocity;
}

bool RigidBodyMotion::moves() const
{
    return !m_position_table.empty() || !m_velocity_table.empty() ||
           !m_orientation_table.empty() || !m_angular_velocity_table.empty() ||
           vs::mag(m_constant_velocity) > 0.0_rt ||
           vs::mag(m_constant_angular_velocity) > 0.0_rt;
}

} // namespace kynema_sgf::actuator::motion
