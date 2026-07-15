#ifndef DRONE_H
#define DRONE_H

#include "src/wind_energy/actuator/actuator_types.H"
#include "src/wind_energy/actuator/sector/ActuatorSector.H"
#include "src/wind_energy/actuator/sector/actuator_sector_ops.H"

#include <memory>
#include <string>

namespace kynema_sgf::actuator {

struct DroneRotor
{
    using SectorData = ActuatorSector::DataType;

    DroneRotor(CFDSim& sim, const std::string& label, int id);

    SectorData data;
    ops::ActSrcOp<ActuatorSector, ActSrcSector> source;
    ops::ProcessOutputsOp<ActuatorSector, ActSrcSector> output;
    vs::Vector body_offset{vs::Vector::zero()};
};

struct DroneData
{
    int num_rotors{0};
    amrex::Real arm_length{0.0_rt};
    RealList arm_lengths;
    RealList arm_angles_degrees;
    amrex::Real arm_phase_degrees{0.0_rt};
    vs::Vector center{vs::Vector::zero()};
    vs::Vector translation_velocity{vs::Vector::zero()};
    vs::Vector body_orientation_degrees{vs::Vector::zero()};
    vs::Tensor body_orientation{vs::Tensor::identity()};
    RealList rotor_omegas;
    amrex::Vector<int> mirror_blades;
    vs::Vector total_force{vs::Vector::zero()};
    vs::Vector total_moment{vs::Vector::zero()};
    amrex::Vector<std::unique_ptr<DroneRotor>> rotors;
};

struct Drone : public ActuatorType
{
    using InfoType = ActInfo;
    using GridType = ActGrid;
    using MetaType = DroneData;
    using DataType = ActDataHolder<Drone>;

    static std::string identifier() { return "Drone"; }
};

namespace drone {

vs::Tensor body_rotation(const vs::Vector& roll_pitch_yaw_degrees);

RealList uniform_arm_angles(int num_rotors, amrex::Real phase_degrees);

VecList
rotor_body_offsets(const RealList& arm_lengths, const RealList& angles_degrees);

} // namespace drone
} // namespace kynema_sgf::actuator

#endif
