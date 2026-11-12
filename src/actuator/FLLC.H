#ifndef FLLC_H
#define FLLC_H

#include <cstdint>
#include <iosfwd>
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
    // width that the simulation resolves where it is above the optimal
    // epsilon: the first value of the input epsilon and, for a turbine, not
    // less than the first value of epsilon_min (read_inputs)
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
 * \param view Force points, velocity points and chord of the line. The two
 * sets of points have to be of one instant, or of instants between which
 * the line moves without changing its shape. The distances along the line
 * taken here are kept for the rest of the run.
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

/**
 * \brief Write what a restart needs of the correction of one line
 *
 * The correction velocity, which the relaxation goes on from, and the
 * distances along the line, which fllc_init took from the positions of the
 * points at the time the correction started (the points of a blade are
 * elsewhere at the restart, and their positions come in single precision).
 * The optimal epsilon and the ratio of the non-uniform points are written
 * to check that the state fits the inputs of the restarted run.
 *
 * \param os Stream of the restart data of the actuator
 * \param data Initialized FLLC data
 */
void fllc_write_state(std::ostream& os, const FLLCData& data);

/**
 * \brief Read the state written by fllc_write_state into initialized data
 *
 * \param is Stream of the restart data of the actuator
 * \param data FLLC data, initialized for the current points and inputs
 * \return False when the state is for other numbers of points, another
 * optimal epsilon or another ratio of the non-uniform points than those of
 * ``data``, or for points whose distances along the line differ from those
 * of ``data`` by more than 2 % of the length of the line. ``data`` is then
 * left as it was
 */
bool fllc_read_state(std::istream& is, FLLCData& data);

} // namespace kynema_sgf::actuator

#endif /* FLLC */
