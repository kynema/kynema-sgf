#include "src/wind_energy/actuator/motion/RotorMotion.H"

#include "src/utilities/trig_ops.H"

#include "AMReX.H"

using namespace amrex::literals;

namespace kynema_sgf::actuator::motion {

void RotorMotion::read_initial_azimuths(
    const utils::ActParser& pp, const int num_rotors)
{
    const bool shared = pp.contains("initial_azimuth_degrees");
    const bool individual = pp.contains("initial_azimuths_degrees");
    if (shared && individual) {
        amrex::Abort(
            "initial_azimuth_degrees and initial_azimuths_degrees are mutually "
            "exclusive");
    }
    m_initial_azimuths.assign(num_rotors, 0.0_rt);
    if (shared) {
        amrex::Real value = 0.0_rt;
        pp.get("initial_azimuth_degrees", value);
        m_initial_azimuths.assign(
            num_rotors, ::kynema_sgf::utils::radians(value));
    } else if (individual) {
        pp.getarr("initial_azimuths_degrees", m_initial_azimuths);
        if (static_cast<int>(m_initial_azimuths.size()) != num_rotors) {
            amrex::Abort(
                "initial_azimuths_degrees must contain num_rotors values");
        }
        for (auto& value : m_initial_azimuths) {
            value = ::kynema_sgf::utils::radians(value);
        }
    }
}

void RotorMotion::read_drone_inputs(
    const utils::ActParser& pp, const int num_rotors)
{
    const bool constant = pp.contains("rotor_omegas");
    const bool timetable = pp.contains("rotor_speed_timetable");
    if (constant == timetable) {
        amrex::Abort(
            "Specify exactly one of rotor_omegas and rotor_speed_timetable");
    }
    m_constant_omegas.assign(num_rotors, 0.0_rt);
    if (constant) {
        pp.getarr("rotor_omegas", m_constant_omegas);
        if (static_cast<int>(m_constant_omegas.size()) != num_rotors) {
            amrex::Abort("rotor_omegas must contain num_rotors values");
        }
    } else {
        std::string filename;
        std::string extrapolation{"hold"};
        pp.get("rotor_speed_timetable", filename);
        pp.query("timetable_extrapolation", extrapolation);
        m_speed_table.read(filename, num_rotors, extrapolation);
    }
    read_initial_azimuths(pp, num_rotors);
}

void RotorMotion::read_sector_inputs(
    const utils::ActParser& pp,
    const amrex::Real preset_omega,
    const bool omega_is_preset)
{
    const bool constant = pp.contains("omega") || omega_is_preset;
    const bool timetable = pp.contains("rotor_speed_timetable");
    if (constant == timetable) {
        amrex::Abort(
            "Specify exactly one of omega and rotor_speed_timetable for an "
            "ActuatorSector");
    }
    m_constant_omegas.assign(1, preset_omega);
    if (!omega_is_preset && pp.contains("omega")) {
        pp.get("omega", m_constant_omegas[0]);
    } else if (timetable) {
        std::string filename;
        std::string extrapolation{"hold"};
        pp.get("rotor_speed_timetable", filename);
        pp.query("timetable_extrapolation", extrapolation);
        m_speed_table.read(filename, 1, extrapolation);
    }
    read_initial_azimuths(pp, 1);
}

amrex::Real RotorMotion::omega(const int rotor, const amrex::Real time) const
{
    AMREX_ALWAYS_ASSERT(rotor >= 0 && rotor < num_rotors());
    return m_speed_table.empty() ? m_constant_omegas[rotor]
                                 : m_speed_table.value(time)[rotor];
}

amrex::Real RotorMotion::azimuth(const int rotor, const amrex::Real time) const
{
    AMREX_ALWAYS_ASSERT(rotor >= 0 && rotor < num_rotors());
    // Integrating omega preserves phase when speed varies; omega(t) * time
    // would only be correct for a constant rotor speed.
    const amrex::Real angle = m_speed_table.empty()
                                  ? m_constant_omegas[rotor] * time
                                  : m_speed_table.integral(time)[rotor];
    return m_initial_azimuths[rotor] + angle;
}

} // namespace kynema_sgf::actuator::motion
