#ifndef ACTUATOR_RIGID_BODY_MOTION_H
#define ACTUATOR_RIGID_BODY_MOTION_H

#include "src/core/vs/quaternion.H"
#include "src/wind_energy/actuator/ActParser.H"
#include "src/wind_energy/actuator/motion/TimeTable.H"

namespace kynema_sgf::actuator::motion {

/** Prescribed translation and orientation of a rigid body.
 *
 * Positions and velocities are in the global CFD frame. The orientation is
 * the active rotation that maps body-frame vectors into that frame.
 */
class RigidBodyMotion
{
public:
    void read_inputs(
        const utils::ActParser& pp,
        const vs::Vector& initial_position,
        const vs::Tensor& initial_orientation,
        const std::string& angular_velocity_key = "angular_velocity",
        const vs::Vector& default_angular_velocity = vs::Vector::zero());

    [[nodiscard]] vs::Vector position(amrex::Real time) const;
    [[nodiscard]] vs::Vector translation_velocity(amrex::Real time) const;
    [[nodiscard]] vs::Tensor orientation(amrex::Real time) const;
    [[nodiscard]] vs::Vector angular_velocity(amrex::Real time) const;

    [[nodiscard]] bool moves() const;

private:
    enum class OrientationFormat { RollPitchYaw, Quaternion };

    [[nodiscard]] vs::Quaternion orientation_quaternion(int row) const;

    vs::Vector m_initial_position{vs::Vector::zero()};
    vs::Vector m_constant_velocity{vs::Vector::zero()};
    vs::Tensor m_initial_orientation{vs::Tensor::identity()};
    vs::Vector m_constant_angular_velocity{vs::Vector::zero()};
    TimeTable m_position_table;
    TimeTable m_velocity_table;
    TimeTable m_orientation_table;
    TimeTable m_angular_velocity_table;
    OrientationFormat m_orientation_format{OrientationFormat::RollPitchYaw};
};

} // namespace kynema_sgf::actuator::motion

#endif
