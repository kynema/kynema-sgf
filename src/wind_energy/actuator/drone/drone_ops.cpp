#include "src/wind_energy/actuator/drone/drone_ops.H"

#include "src/CFDSim.H"
#include "src/utilities/constants.H"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>

#include "AMReX_Print.H"

using namespace amrex::literals;

namespace kynema_sgf::actuator::ops {
namespace {

void input_error(const std::string& label, const std::string& message)
{
    amrex::Abort("Drone '" + label + "': " + message);
}

void validate_distinct_angles(const std::string& label, const RealList& angles)
{
    for (int i = 0; i < static_cast<int>(angles.size()); ++i) {
        for (int j = i + 1; j < static_cast<int>(angles.size()); ++j) {
            amrex::Real delta = std::remainder(angles[i] - angles[j], 360.0_rt);
            if (std::abs(delta) < constants::TIGHT_TOL) {
                input_error(label, "arm angles must be distinct");
            }
        }
    }
}

} // namespace

void ReadInputsOp<Drone, ActSrcDrone>::operator()(
    Drone::DataType& data, const utils::ActParser& pp)
{
    auto& meta = data.meta();
    const auto& label = data.info().label;

    pp.get("num_rotors", meta.num_rotors);
    if (meta.num_rotors <= 0) {
        input_error(label, "num_rotors must be positive");
    }
    if (!pp.contains("position_timetable")) {
        pp.get("center", meta.center);
    }
    pp.query("translation_velocity", meta.translation_velocity);
    pp.query("body_orientation_degrees", meta.body_orientation_degrees);
    if (pp.contains("orientation_timetable") &&
        pp.contains("body_orientation_degrees")) {
        input_error(
            label,
            "body_orientation_degrees cannot be combined with "
            "orientation_timetable");
    }
    meta.body_orientation = drone::body_rotation(meta.body_orientation_degrees);
    // All child sectors share these histories, keeping the rigid-body pose and
    // rotor clocks synchronized across the composite actuator.
    meta.body_motion = std::make_shared<motion::RigidBodyMotion>();
    meta.body_motion->read_inputs(
        pp, meta.center, meta.body_orientation, "angular_velocity");
    meta.rotor_motion = std::make_shared<motion::RotorMotion>();
    meta.rotor_motion->read_drone_inputs(pp, meta.num_rotors);

    // One arm length applies to every rotor; num_rotors values permit an
    // asymmetric layout. MultiParser preserves instance-over-default priority.
    if (!pp.contains("arm_length")) {
        input_error(label, "arm_length is required");
    }
    pp.getarr("arm_length", meta.arm_lengths);
    if (meta.arm_lengths.size() == 1) {
        const auto arm_length = meta.arm_lengths.front();
        meta.arm_lengths.assign(meta.num_rotors, arm_length);
    }
    if (static_cast<int>(meta.arm_lengths.size()) != meta.num_rotors) {
        input_error(
            label,
            "arm_length must contain either one value or num_rotors values");
    }
    if (std::ranges::any_of(meta.arm_lengths, [](const amrex::Real length) {
            return length <= 0.0_rt;
        })) {
        input_error(label, "all arm lengths must be positive");
    }

    // Explicit angles describe irregular layouts. Otherwise phase rotates a
    // uniformly spaced layout about the body z axis.
    const auto& instance = pp.params();
    const auto& defaults = pp.default_params();
    const bool instance_angles = instance.contains("arm_angles_degrees");
    const bool instance_phase = instance.contains("arm_phase_degrees");
    const bool instance_layout = instance_angles || instance_phase;
    const bool conflicting_layout =
        (instance_angles && instance_phase) ||
        (!instance_layout && defaults.contains("arm_angles_degrees") &&
         defaults.contains("arm_phase_degrees"));
    if (conflicting_layout) {
        input_error(
            label,
            "arm_angles_degrees and arm_phase_degrees are mutually exclusive");
    }
    if (instance_angles) {
        pp.getarr("arm_angles_degrees", meta.arm_angles_degrees);
    } else if (instance_phase) {
        pp.get("arm_phase_degrees", meta.arm_phase_degrees);
        meta.arm_angles_degrees =
            drone::uniform_arm_angles(meta.num_rotors, meta.arm_phase_degrees);
    } else if (defaults.contains("arm_angles_degrees")) {
        pp.getarr("arm_angles_degrees", meta.arm_angles_degrees);
    } else {
        pp.query("arm_phase_degrees", meta.arm_phase_degrees);
        meta.arm_angles_degrees =
            drone::uniform_arm_angles(meta.num_rotors, meta.arm_phase_degrees);
    }
    if (static_cast<int>(meta.arm_angles_degrees.size()) != meta.num_rotors) {
        input_error(label, "arm_angles_degrees must contain num_rotors values");
    }
    validate_distinct_angles(label, meta.arm_angles_degrees);

    meta.mirror_blades.assign(meta.num_rotors, false);
    if (pp.contains("mirror_blades")) {
        amrex::Vector<std::string> mirror_inputs;
        pp.getarr("mirror_blades", mirror_inputs);
        if (static_cast<int>(mirror_inputs.size()) != meta.num_rotors) {
            input_error(label, "mirror_blades must contain num_rotors values");
        }
        for (int i = 0; i < meta.num_rotors; ++i) {
            const auto value = amrex::toLower(mirror_inputs[i]);
            if ((value != "true") && (value != "false")) {
                input_error(
                    label, "mirror_blades values must be true or false");
            }
            meta.mirror_blades[i] = (value == "true");
        }
    }

    const auto offsets =
        drone::rotor_body_offsets(meta.arm_lengths, meta.arm_angles_degrees);
    meta.rotors.reserve(meta.num_rotors);
    for (int i = 0; i < meta.num_rotors; ++i) {
        const std::string rotor_label = label + ".R" + std::to_string(i + 1);
        auto rotor = std::make_unique<DroneRotor>(
            data.sim(), rotor_label, data.info().id * 100000 + i);
        rotor->body_offset = offsets[i];

        // A drone owns its rotor configuration. Reuse the same default and
        // instance namespaces for both drone geometry and sector aerodynamics;
        // each reader simply ignores inputs outside its responsibility.
        utils::ActParser sector_pp("Actuator.Drone", "Actuator." + label);
        auto& rotor_meta = rotor->data.meta();
        rotor_meta.body_motion = meta.body_motion;
        rotor_meta.rotor_motion = meta.rotor_motion;
        rotor_meta.rotor_index = i;
        rotor_meta.body_offset = offsets[i];
        rotor_meta.omega = meta.rotor_motion->omega(i, 0.0_rt);
        rotor_meta.user_omega = true;
        ReadInputsOp<ActuatorSector, ActSrcSector>()(rotor->data, sector_pp);
        if (meta.mirror_blades[i]) {
            // Mirroring the blade geometry reverses twist without introducing
            // a separate rotation-direction convention.
            for (auto& twist : rotor->data.meta().twist_inp) {
                twist = -twist;
            }
        }
        const auto rotor_center =
            meta.body_motion->position(0.0_rt) +
            (meta.body_motion->orientation(0.0_rt) & rotor->body_offset);
        const auto initial_normal =
            meta.body_motion->orientation(0.0_rt) & vs::Vector::khat();
        sector::set_placement(
            rotor->data, rotor_center, initial_normal,
            meta.body_motion->translation_velocity(0.0_rt));
        rotor->output.read_io_options(sector_pp);
        meta.rotors.emplace_back(std::move(rotor));
    }

    // The composite search region is the union of its child-sector regions.
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> lo;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> hi;
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        lo[n] = std::numeric_limits<amrex::Real>::max();
        hi[n] = std::numeric_limits<amrex::Real>::lowest();
    }
    for (const auto& rotor : meta.rotors) {
        const auto& box = rotor->data.info().bound_box;
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            lo[n] = std::min(lo[n], box.lo(n));
            hi[n] = std::max(hi[n], box.hi(n));
        }
    }
    data.info().bound_box = amrex::RealBox(lo.data(), hi.data());
}

void InitDataOp<Drone, ActSrcDrone>::operator()(Drone::DataType& data)
{
    auto& meta = data.meta();
    int total_points = 0;
    for (auto& rotor : meta.rotors) {
        InitDataOp<ActuatorSector, ActSrcSector>()(rotor->data);
        total_points += static_cast<int>(rotor->data.grid().vel_pos.size());
    }
    data.grid().resize(0, total_points);
    UpdatePosOp<Drone, ActSrcDrone>()(data);
}

void UpdatePosOp<Drone, ActSrcDrone>::operator()(Drone::DataType& data)
{
    auto& grid = data.grid();
    int offset = 0;
    for (auto& rotor : data.meta().rotors) {
        UpdatePosOp<ActuatorSector, ActSrcSector>()(rotor->data);
        const auto& positions = rotor->data.grid().vel_pos;
        std::copy(
            positions.begin(), positions.end(), grid.vel_pos.begin() + offset);
        offset += static_cast<int>(positions.size());
    }
}

void UpdateVelOp<Drone, ActSrcDrone>::operator()(Drone::DataType& data)
{
    auto& grid = data.grid();
    int offset = 0;
    for (auto& rotor : data.meta().rotors) {
        auto& rotor_grid = rotor->data.grid();
        const int npts = static_cast<int>(rotor_grid.vel.size());
        std::copy_n(grid.vel.begin() + offset, npts, rotor_grid.vel.begin());
        std::copy_n(
            grid.density.begin() + offset, npts, rotor_grid.density.begin());
        offset += npts;
    }
}

void ComputeForceOp<Drone, ActSrcDrone>::operator()(Drone::DataType& data)
{
    auto& meta = data.meta();
    const auto& time = data.sim().time();
    const amrex::Real midpoint_time =
        time.current_time() + 0.5_rt * sector::timestep_width(data.sim());
    const auto drone_center = meta.body_motion->position(midpoint_time);
    meta.total_force = vs::Vector::zero();
    meta.total_moment = vs::Vector::zero();
    for (auto& rotor : meta.rotors) {
        ComputeForceOp<ActuatorSector, ActSrcSector>()(rotor->data);
        const auto& rotor_meta = rotor->data.meta();
        // Sector forces act on the fluid; report equal-and-opposite vehicle
        // load.
        meta.total_force = meta.total_force - rotor_meta.integrated_force;
        meta.total_moment =
            meta.total_moment - rotor_meta.integrated_moment -
            ((rotor_meta.center - drone_center) ^ rotor_meta.integrated_force);
    }
}

void ProcessOutputsOp<Drone, ActSrcDrone>::prepare_outputs(
    const std::string& out_dir)
{
    m_filename = out_dir + "/" + m_data.info().label + "_loads.csv";
    std::ofstream stream(m_filename);
    stream << "time,force_x,force_y,force_z,moment_x,moment_y,moment_z\n";
    for (auto& rotor : m_data.meta().rotors) {
        rotor->output.prepare_outputs(out_dir);
    }
}

void ProcessOutputsOp<Drone, ActSrcDrone>::write_outputs()
{
    const auto& time = m_data.sim().time();
    if ((m_output_frequency > 0) &&
        (time.time_index() % m_output_frequency == 0)) {
        const auto& force = m_data.meta().total_force;
        const auto& moment = m_data.meta().total_moment;
        std::ofstream stream(m_filename, std::ios::app);
        stream << time.new_time() << ',' << force.x() << ',' << force.y() << ','
               << force.z() << ',' << moment.x() << ',' << moment.y() << ','
               << moment.z() << '\n';
    }
    for (auto& rotor : m_data.meta().rotors) {
        rotor->output.write_outputs();
    }
}

void ActSrcOp<Drone, ActSrcDrone>::initialize()
{
    for (auto& rotor : m_data.meta().rotors) {
        rotor->source.initialize();
    }
}

void ActSrcOp<Drone, ActSrcDrone>::setup_op()
{
    for (auto& rotor : m_data.meta().rotors) {
        rotor->source.setup_op();
    }
}

void ActSrcOp<Drone, ActSrcDrone>::operator()(
    const int lev, const amrex::MFIter& mfi, const amrex::Geometry& geom)
{
    for (auto& rotor : m_data.meta().rotors) {
        rotor->source(lev, mfi, geom);
    }
}

} // namespace kynema_sgf::actuator::ops
