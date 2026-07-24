#ifndef ACTUATORSECTOR_H
#define ACTUATORSECTOR_H

#include "src/wind_energy/actuator/actuator_types.H"
#include "src/wind_energy/actuator/aero/AirfoilTable.H"
#include "src/wind_energy/actuator/motion/RigidBodyMotion.H"
#include "src/wind_energy/actuator/motion/RotorMotion.H"

#include <memory>
#include <string>

#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::actuator {

/** Persistent data for the actuator-sector rotor model.
 *
 *  Frame convention:
 *  - CFD/fixed frame: all vectors stored in ActGrid and motion vectors.
 *  - Rotor local frame: x is radial at zero azimuth, y is clockwise-positive
 *    tangential direction at zero azimuth, and z is the rotor normal.
 *  - Positive omega advances a blade clockwise when viewed along local +z,
 *    matching the theory/prototype convention:
 *      e_r = (cos(theta), -sin(theta), 0)
 *      e_theta = (-sin(theta), -cos(theta), 0)
 */
struct ActuatorSectorData
{
    //! Prescribed rigid-body motion shared with an owning Drone when present
    std::shared_ptr<motion::RigidBodyMotion> body_motion;

    //! Prescribed rotor-speed motion shared with an owning Drone when present
    std::shared_ptr<motion::RotorMotion> rotor_motion;

    //! Rotor index within rotor_motion
    int rotor_index{0};

    //! Hub offset from the rigid-body center in body coordinates [m]
    vs::Vector body_offset{vs::Vector::zero()};

    //! Number of rotor blades [-]
    int num_blades{2};

    //! Rotor diameter [m]
    amrex::Real rotor_diameter{1.0_rt};

    //! Rotor radius [m]
    amrex::Real rotor_radius{0.5_rt};

    //! Inboard cutoff as a fraction of rotor radius [-]
    amrex::Real root_radius_fraction{0.18_rt};

    //! Inboard cutoff radius [m]
    amrex::Real root_radius{0.09_rt};

    //! Rotor spin angular speed [rad/s]
    amrex::Real omega{0.0_rt};
    bool user_omega{false};

    //! Initial rotor hub center in the CFD frame [m]
    vs::Vector center0{0.0_rt, 0.0_rt, 0.0_rt};

    //! Rotor hub center at the current force-evaluation midpoint [m]
    vs::Vector center{0.0_rt, 0.0_rt, 0.0_rt};

    //! Constant drone/body translation velocity in the CFD frame [m/s]
    vs::Vector translation_velocity{0.0_rt, 0.0_rt, 0.0_rt};

    //! Initial rotor normal direction in the CFD frame [-]
    vs::Vector rotor_normal{1.0_rt, 0.0_rt, 0.0_rt};

    //! Rotor-frame rotation rate specified per blade revolution [deg/rev]
    vs::Vector rotor_rotation_degrees_per_revolution{0.0_rt, 0.0_rt, 0.0_rt};

    //! Rotor-frame angular velocity in the CFD frame [rad/s]
    vs::Vector rotor_angular_velocity{0.0_rt, 0.0_rt, 0.0_rt};

    //! True when rotor_angular_velocity was explicitly supplied by the user
    bool user_rotor_angular_velocity{false};

    //! Rotor azimuth at the current force-evaluation midpoint [rad]
    amrex::Real azimuth{0.0_rt};

    //! Blade azimuth swept over the current timestep [rad]
    amrex::Real delta_azimuth{0.0_rt};

    //! Non-dimensional locations for chord/twist interpolation [-]
    RealList span_locs{0.0_rt, 1.0_rt};

    //! Input chord distribution at span_locs [m]
    RealList chord_inp{1.0_rt, 1.0_rt};

    //! Input twist distribution at span_locs [deg]
    RealList twist_inp{0.0_rt, 0.0_rt};

    //! Airfoil table filename
    std::string airfoil_file;

    //! Airfoil table format, e.g., openfast
    std::string airfoil_type{"openfast"};

    //! Sectional airfoil lookup table, returning Cl/Cd as a function of AoA
    std::unique_ptr<AirfoilTable> aflookup;

    //! Fixed Gaussian width [m]
    amrex::Real epsilon{0.0_rt};

    //! True when epsilon was explicitly supplied by the user
    bool user_epsilon{false};

    //! Gaussian width divided by local chord [-]
    amrex::Real epsilon_chord{0.0_rt};

    //! Minimum Gaussian width when using epsilon_chord [m]
    amrex::Real epsilon_min{0.0_rt};

    //! Required local span resolution, epsilon / dr >= epsilon_dr [-]
    amrex::Real epsilon_dr{1.5_rt};

    //! Required swept arc length resolution, epsilon / dl >= epsilon_dl [-]
    amrex::Real epsilon_dl{1.5_rt};

    //! Required local chord resolution, chord / dr >= min_chord_dr [-]
    amrex::Real min_chord_dr{2.0_rt};

    //! Compact Gaussian support radius divided by epsilon [-]
    amrex::Real support_radius_over_epsilon{3.0_rt};

    //! Span quadrature rule. Initial implementation supports midpoint.
    std::string radial_quadrature{"midpoint"};

    //! Swept quadrature rule. Initial implementation supports midpoint.
    std::string swept_quadrature{"midpoint"};

    //! Gaussian evaluation mode: direct or table
    std::string gaussian_type{"direct"};

    //! Maximum interpolation error for the Gaussian exp(-r^2) table [-]
    amrex::Real gaussian_table_error{1.0e-4_rt};

    //! Number of interpolation intervals in the Gaussian table [-]
    int gaussian_table_nintervals{0};

    //! Host-side table for exp(-r^2), r normalized by epsilon [-]
    RealList gaussian_table;

    //! Radial blade station locations [m]
    RealList radius;

    //! Radial quadrature weights [m]
    RealList dr;

    //! Local chord at radial stations [m]
    RealList chord;

    //! Local twist at radial stations [deg]
    RealList twist;

    //! Local Gaussian width at radial stations [m]
    RealList epsilon_profile;

    //! Number of swept quadrature samples at each radial station [-]
    amrex::Vector<int> theta_counts;

    //! Relative velocity in the local blade-section plane [m/s]
    VecList vel_rel;

    //! Force on fluid at midpoint blade stations in fixed-frame coordinates
    //! before multiplication by radial quadrature weight [m^3/s^2]
    VecList blade_force;

    //! Force on fluid along local tangential direction at blade stations
    //! before multiplication by radial quadrature weight [m^3/s^2]
    RealList blade_force_theta;

    //! Force on fluid along local rotor normal at blade stations before
    //! multiplication by radial quadrature weight [m^3/s^2]
    RealList blade_force_normal;

    //! Angle of attack at midpoint blade stations [deg]
    RealList aoa;

    //! Lift coefficient at midpoint blade stations [-]
    RealList cl;

    //! Drag coefficient at midpoint blade stations [-]
    RealList cd;

    //! Integrated force component along the rotor normal [m^4/s^2]
    amrex::Real thrust{0.0_rt};

    //! Integrated torque about the rotor normal [m^5/s^2]
    amrex::Real torque{0.0_rt};

    //! Integrated density-normalized force applied to the fluid in the CFD
    //! frame [m^4/s^2]
    vs::Vector integrated_force{0.0_rt, 0.0_rt, 0.0_rt};

    //! Integrated density-normalized moment applied to the fluid about the
    //! rotor hub [m^5/s^2]
    vs::Vector integrated_moment{0.0_rt, 0.0_rt, 0.0_rt};
};

struct ActuatorSector : public ActuatorType
{
    using InfoType = ActInfo;
    using GridType = ActGrid;
    using MetaType = ActuatorSectorData;
    using DataType = ActDataHolder<ActuatorSector>;

    static std::string identifier() { return "Actuator"; }
};

} // namespace kynema_sgf::actuator

#endif /* ACTUATORSECTOR_H */
