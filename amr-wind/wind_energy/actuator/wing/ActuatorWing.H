#ifndef ACTUATORWING_H
#define ACTUATORWING_H

#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/core/vs/vector_space.H"

#include <string>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator {

struct FLLCData;

/** Base data representation of a fixed-wing.
 *
 *  \ingroup actuator
 */
struct WingBaseData
{
    //! Number of points along the wing
    int num_pts;

    //! Starting coordinate of the wing
    vs::Vector start;

    //! Ending coordinate of the wing
    vs::Vector end;

    //! The normal vector perpendicular to the span
    vs::Vector blade_x{1.0_rt, 0, 0};

    //! Gaussian smearing factor input by user
    vs::Vector eps_inp;

    //! epsilon / chord
    vs::Vector epsilon_chord;

    //! Pitch angle for the wing
    amrex::Real pitch{0.0_rt};

    //! Switch for turning on 2D Gaussian (off in spanwise)
    bool gauss_2D{false};

    //! Switch for spanwise normalization in 2D setting
    bool normalize_2D_spanwise{true};

    //! Switch for turning off force components of actuator
    amrex::Vector<int> force_coord_flags{1, 1, 1};

    //! User-specified velocity magnitude, only active when nonnegative
    amrex::Real prescribed_uinf{-1.0_rt};

    //! File name for pitch actuation table
    std::string pitch_timetable_file;

    //! Pitch actuation: time table
    amrex::Vector<amrex::Real> time_table;

    //! Pitch actuation: pitch angle table
    amrex::Vector<amrex::Real> pitch_table;

    //! Relative velocity at the actuator node
    VecList vel_rel;

    //! Translational velocity
    vs::Vector vel_tr{0.0_rt, 0.0_rt, 0.0_rt};

    //! Wing motion type
    std::string motion_type{"none"};

    //! Parameters for sine motions
    vs::Vector s_vector{0.0_rt, 0.0_rt, 0.0_rt};
    amrex::Real s_period{0.0_rt};

    //! Filtered Lifting Line Correction data
    std::unique_ptr<FLLCData> fllc;

    //! Data structure required for passing data to the FLLC
    ComponentView component_view;

    //! Width of actuator sections
    RealList dx;

    //! Angle of attack at actuator nodes
    RealList aoa;

    //! Computed lift coefficient at actuator nodes
    RealList cl;

    //! Computed drag coefficient at actuator nodes
    RealList cd;

    //! Chord length at the actuator nodes
    RealList chord;

    //! Total integrated lift over the wing
    amrex::Real lift{0.0_rt};

    //! Total integrated drag over the wing
    amrex::Real drag{0.0_rt};
};

struct WingType : public ActuatorType
{};

} // namespace amr_wind::actuator

#endif /* ACTUATORWING_H */
