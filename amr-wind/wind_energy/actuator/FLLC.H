#ifndef FLLC_H
#define FLLC_H

#include <cstdint>
#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator {

enum class FLLCType : std::uint8_t { ConstantChord, VariableChord };

static std::map<std::string, FLLCType> FLLCTypeMap{
    {"constant_chord", FLLCType::ConstantChord},
    {"variable_chord", FLLCType::VariableChord}};

struct FLLCData
{
    // constants
    FLLCType correction_type{FLLCType::VariableChord};
    amrex::Real epsilon;
    amrex::Real relaxation_factor{0.1_rt};
    amrex::Real fllc_start_time{0.0_rt};
    RealList r;
    RealList dr;
    RealList optimal_epsilon;

    // computed values
    VecList force_point_velocity;
    VecList les_velocity;
    VecList optimal_velocity;
    VecList correction_velocity;
    VecList lift;
    VecList grad_lift;
    bool initialized{false};
    bool different_sizes;
    RealList span_distance_vel;
    RealList span_distance_force;

    // non-uniform variables
    bool nonuniform{true};      // non-uniform flag
    amrex::Real eps_dr{1.};     // the ratio of epsilon to actuator width
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
 * \param pp Parser
 * \param data Data to be populated from parser
 */
void fllc_parse(const utils::ActParser& pp, FLLCData& data);

/**
 * \brief Initialize FLLC data structure. This should be called at the end of
 * the first ComputeForceOp to ensure the data is fully populated
 *
 * \param data FLLC specific data and scratch fields
 * \param view Component view that has references to the grid data
 * \param eps_chord epsilon chord
 */
void fllc_init(
    FLLCData& data, const ComponentView& view, amrex::Real eps_chord);

} // namespace amr_wind::actuator

#endif /* FLLC */
