#ifndef EXTERNAL_BASE_TYPES_H
#define EXTERNAL_BASE_TYPES_H

#include <string>
#include <array>
#include <cstdint>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace ext_turb {

enum class SimMode : std::uint8_t {
    init = 0, // Clean start
    replay,   // Replay using velocities stored in file
    restart   // Restart using external checkpoint files
};

template <typename ExtData>
std::string ext_id();

/* Local axes of the blade points of an external solver: rows of the
 *  orientation matrix it sends (grid.orientation after scatter_data)
 */
enum class BladeFrame : std::uint8_t {
    // x normal to chord, y along chord, z along span (OpenFAST: AeroDyn
    // airfoil frame)
    normal_chord_span = 0,
    // x along span, y along chord, z normal to chord (Kynema-FMB: aero
    // sections)
    span_chord_normal
};

/* Blade frame of a solver turbine type: its static member blade_frame,
 *  else the OpenFAST frame (test stand-ins)
 */
template <typename SolverTurbine>
constexpr BladeFrame blade_frame_of()
{
    if constexpr (requires { SolverTurbine::blade_frame; }) {
        return SolverTurbine::blade_frame;
    } else {
        return BladeFrame::normal_chord_span;
    }
}

/* Stored epsilon component of input value n
 *
 *  Inputs: (chord, thickness, span). Stored: local frame of the point,
 *  component m along row m of its orientation.
 */
constexpr int epsilon_component(const BladeFrame frame, const int n)
{
    switch (frame) {
    case BladeFrame::span_chord_normal:
        return (n == 0) ? 1 : ((n == 1) ? 2 : 0);
    case BladeFrame::normal_chord_span:
    default:
        return (n == 0) ? 1 : ((n == 1) ? 0 : 2);
    }
}

/* Representation of a turbine for exchanging data with external solver */
struct ExternalTurbine
{
    static constexpr int ndim = 3;

    // Unique string identifier for this turbine
    std::string tlabel;

    // Local ID for this turbine (provided by external interface)
    int tid_local;

    // Global ID for this turbine (set by Turbine instance)
    int tid_global;

    // Number of actuator points per blade
    int num_pts_blade;

    // Number of actuator points for tower
    int num_pts_tower;

    // Position of tower base in global coordinate system
    std::array<float, ndim> base_pos;

    // Hub stats for Spinner Lidar
    std::array<float, 3> hub_abs_pos;

    std::array<float, 3> hub_rot_vel;

    std::array<amrex::Real, 9> hub_orient;

    SimMode sim_mode{SimMode::init};

    // External input file
    std::string input_file;

    // Checkpoint file name
    std::string checkpoint_file;

    // Number of blades
    int num_blades;

    // Total number of elements along the blade
    int num_blade_elem;
    // Total number of elements along the tower
    int num_tower_elem;

    bool is_solution0{true};

    // Start time for this turbine
    amrex::Real start_time{0.0_rt};

    // End time for this turbine
    amrex::Real stop_time{-1.0_rt};

    // Timestep for CFD
    amrex::Real dt_cfd;

    // Timestep for external turbine solver
    amrex::Real dt_ext;

    // Number of sub-steps of external solver per CFD timestep
    int num_substeps;

    // Time step index for external
    int time_index{0};

    // Data access functions that have to be defined for each type
    [[nodiscard]] virtual float* position_at_vel(int dir) const = 0;
    [[nodiscard]] virtual float* solid_velocity(int dir) const = 0;
    [[nodiscard]] virtual float* fluid_velocity(int dir) const = 0;
    [[nodiscard]] virtual float* force(int dir) const = 0;
    [[nodiscard]] virtual float* position_at_force(int dir) const = 0;
    [[nodiscard]] virtual float* chord_at_force() const = 0;
    [[nodiscard]] virtual float* orientation() const = 0;

    [[nodiscard]] virtual int length_fluid_velocity(int dir) const = 0;
    [[nodiscard]] virtual int length_force(int dir) const = 0;
    [[nodiscard]] virtual int length_position_at_force(int dir) const = 0;
    [[nodiscard]] virtual int length_orientation() const = 0;
    [[nodiscard]] virtual int num_vel_pts_blade() const = 0;
};

} // namespace ext_turb

#endif /* EXTERNAL_BASE_TYPES_H */
