#include "src/wind_energy/actuator/drone/Drone.H"
#include "src/wind_energy/actuator/drone/drone_ops.H"
#include "src/wind_energy/actuator/ActuatorModel.H"
#include "src/utilities/trig_ops.H"

#include <cmath>

namespace kynema_sgf::actuator {

DroneRotor::DroneRotor(CFDSim& sim, const std::string& label, const int id)
    : data(sim, label, id), source(data), output(data)
{}

namespace drone {

vs::Tensor body_rotation(const vs::Vector& angles)
{
    // Apply intrinsic body x, then y, then z rotations to body-frame vectors.
    return vs::zrot(angles.z()) & vs::yrot(angles.y()) & vs::xrot(angles.x());
}

RealList uniform_arm_angles(const int num_rotors, const amrex::Real phase)
{
    RealList angles(num_rotors);
    for (int i = 0; i < num_rotors; ++i) {
        angles[i] = phase + 360.0_rt * static_cast<amrex::Real>(i) /
                                static_cast<amrex::Real>(num_rotors);
    }
    return angles;
}

VecList
rotor_body_offsets(const RealList& lengths, const RealList& angles_degrees)
{
    AMREX_ALWAYS_ASSERT(lengths.size() == angles_degrees.size());
    // Arms lie in the body x-y plane; rigid-body motion maps these fixed
    // offsets into the CFD frame at runtime.
    VecList offsets(lengths.size());
    for (int i = 0; i < static_cast<int>(lengths.size()); ++i) {
        const amrex::Real angle =
            ::kynema_sgf::utils::radians(angles_degrees[i]);
        offsets[i] = {
            lengths[i] * std::cos(angle), lengths[i] * std::sin(angle), 0.0_rt};
    }
    return offsets;
}

} // namespace drone

template class ActModel<Drone, ActSrcDrone>;

} // namespace kynema_sgf::actuator
