#ifndef EXTERNAL_BASE_TYPES_H
#define EXTERNAL_BASE_TYPES_H

#include <string>
#include <cstdint>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace ext_turb {

enum class SimMode : std::uint8_t {
    init = 0, ///< Clean start
    replay,   ///< Replay using velocities stored in file
    restart   ///< Restart using external checkpoint files
};

template <typename ExtData>
std::string ext_id();

/** Representation of a turbine for exchanging data with external solver
 */
struct ExternalTurbine
{
    static constexpr int ndim = 3;

    //! Unique string identifier for this turbine
    std::string tlabel;

    //! Local ID for this turbine (provided by external interface)
    int tid_local;

    //! Global ID for this turbine (set by Turbine instance)
    int tid_global;

    //! Number of actuator points per blade
    int num_pts_blade;

    //! Number of actuator points for tower
    int num_pts_tower;

    //! Position of tower base in global coordinate system
    amrex::Array<float, ndim> base_pos;

    //! Hub stats for Spinner Lidar
    amrex::Array<float, 3> hub_abs_pos;

    amrex::Array<float, 3> hub_rot_vel;

    amrex::Array<amrex::Real, 9> hub_orient;

    SimMode sim_mode{SimMode::init};

    //! External input file
    std::string input_file;

    //! Checkpoint file name
    std::string checkpoint_file;

    //! Number of blades
    int num_blades;

    //! Total number of elements along the blade
    int num_blade_elem;
    //! Total number of elements along the tower
    int num_tower_elem;

    bool is_solution0{true};

    //! Start time for this turbine
    amrex::Real start_time{0.0_rt};

    //! End time for this turbine
    amrex::Real stop_time{-1.0_rt};

    //! Timestep for CFD
    amrex::Real dt_cfd;

    //! Timestep for external turbine solver
    amrex::Real dt_ext;

    //! Number of sub-steps of external solver per CFD timestep
    int num_substeps;

    //! Time step index for external
    int time_index{0};

    //! Checkpoint interval for external
    int chkpt_interval;

    //! Data access functions that have to be defined for each type
    virtual float* position_at_vel(int dir) const = 0;
    virtual float* solid_velocity(int dir) const = 0;
    virtual float* fluid_velocity(int dir) const = 0;
    virtual float* force(int dir) const = 0;
    virtual float* position_at_force(int dir) const = 0;
    virtual float* chord_at_force() const = 0;
    virtual float* orientation() const = 0;

    virtual int length_fluid_velocity(int dir) const = 0;
    virtual int length_force(int dir) const = 0;
    virtual int length_position_at_force(int dir) const = 0;
    virtual int length_orientation() const = 0;
    virtual int num_vel_pts_blade() const = 0;
};

} // namespace ext_turb

#endif /* EXTERNAL_BASE_TYPES_H */
