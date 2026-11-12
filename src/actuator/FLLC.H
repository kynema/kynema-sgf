#ifndef FLLC_H
#define FLLC_H

#include <cstdint>
#include <string>
#include "src/actuator/actuator_types.H"
#include "src/actuator/ActParser.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::actuator {

enum class FLLCType : std::uint8_t { ConstantChord, VariableChord };

static std::map<std::string, FLLCType> FLLCTypeMap{
    {"constant_chord", FLLCType::ConstantChord},
    {"variable_chord", FLLCType::VariableChord}};

struct FLLCData
{
    // constants
    FLLCType correction_type{FLLCType::VariableChord};
    amrex::Real epsilon{0.0_rt};
    amrex::Real relaxation_factor{0.1_rt};
    amrex::Real fllc_start_time{0.0_rt};
    RealList r;
    RealList dr;
    RealList optimal_epsilon;

    // computed values
    VecList les_velocity;
    VecList optimal_velocity;
    VecList correction_velocity;
    VecList lift;
    VecList grad_lift;
    bool initialized{false};
    // true when the velocity points are not at the force points
    bool mismatched_points{false};
    // distance of each velocity/force point from the first one of its set
    RealList span_distance_vel;
    RealList span_distance_force;

    // non-uniform variables
    bool nonuniform{true};      // non-uniform flag
    amrex::Real eps_dr{3.0_rt}; // the ratio of epsilon to actuator width
    RealList nonuniform_r;      // non-uniform radius
    RealList nonuniform_dr;     // non-uniform spacing
    VecList vel_rel;            // uniform relative velocity
    VecList nonuniform_vel_rel; // non-uniform relative velocity
    RealList nonuniform_optimal_epsilon; // non-uniform radius
    VecList nonuniform_lift;             // non-uniform lift
};

/**
 * \brief Function to capture common parsing requirements for the filtered
 * lifting line correction
 *
 * Aborts, naming the input, when fllc_relaxation_factor is not in (0, 1] or
 * fllc_epsilon_dr_ratio is not larger than 0.
 *
 * \param pp Parser
 * \param data Data to be populated from parser
 */
void fllc_parse(const utils::ActParser& pp, FLLCData& data);

/**
 * \brief Initialize FLLC data structure. This should be called at the end of
 * the first ComputeLoadsOp to ensure the data is fully populated
 *
 * Aborts when the line cannot be used for the correction (fewer than two
 * points, points that do not move away from the first one, a chord or an
 * epsilon_chord that is not larger than 0). Warns when the actuator points
 * do not resolve the optimal epsilon.
 *
 * \param data FLLC specific data and scratch fields
 * \param view Component view that has references to the grid data
 * \param eps_chord epsilon chord
 * \param label Name of the line (actuator, blade) for the messages
 * \param print_warnings False on the processes that are to stay silent
 * (every process of the actuator initializes the correction)
 */
void fllc_init(
    FLLCData& data,
    const ComponentView& view,
    amrex::Real eps_chord,
    const std::string& label = "",
    bool print_warnings = true);

} // namespace kynema_sgf::actuator

#endif /* FLLC */
