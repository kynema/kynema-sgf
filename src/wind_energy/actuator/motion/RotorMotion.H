#ifndef ACTUATOR_ROTOR_MOTION_H
#define ACTUATOR_ROTOR_MOTION_H

#include "src/wind_energy/actuator/ActParser.H"
#include "src/wind_energy/actuator/motion/TimeTable.H"

namespace kynema_sgf::actuator::motion {

/** Signed rotor-speed histories and their corresponding azimuths. */
class RotorMotion
{
public:
    void read_drone_inputs(const utils::ActParser& pp, int num_rotors);
    void read_sector_inputs(
        const utils::ActParser& pp,
        amrex::Real preset_omega,
        bool omega_is_preset);

    [[nodiscard]] int num_rotors() const
    {
        return static_cast<int>(m_constant_omegas.size());
    }
    [[nodiscard]] amrex::Real omega(int rotor, amrex::Real time) const;
    [[nodiscard]] amrex::Real azimuth(int rotor, amrex::Real time) const;

private:
    void read_initial_azimuths(const utils::ActParser& pp, int num_rotors);

    RealList m_constant_omegas;
    //! Initial rotor angles stored internally in radians
    RealList m_initial_azimuths;
    TimeTable m_speed_table;
};

} // namespace kynema_sgf::actuator::motion

#endif
