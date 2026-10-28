#ifndef KYNEMA_FMB_HOST_H
#define KYNEMA_FMB_HOST_H

#include <memory>
#include <string>
#include "AMReX_REAL.H"

/** Boundary between Kynema-SGF and the Kynema-FMB library
 *
 *  nvcc does not lay out Kokkos::View the same way as the host compiler, so
 *  all code that uses Kynema-FMB types lives in KynemaFMBHost.cpp, which is
 *  compiled with the same compiler as Kokkos (see CMakeLists.txt). With a
 *  host-only Kokkos that file cannot include AMReX headers containing GPU
 *  code, so this header only uses plain C++ types and Kynema-FMB types are
 *  only forward declared.
 */
namespace kynema_fmb::interfaces {
class TurbineInterface;
} // namespace kynema_fmb::interfaces

namespace sgf_fmb {

using TurbineInterface = kynema_fmb::interfaces::TurbineInterface;

struct TurbineInterfaceDeleter
{
    void operator()(TurbineInterface* ptr) const;
};

using TurbineInterfacePtr =
    std::unique_ptr<TurbineInterface, TurbineInterfaceDeleter>;

//! Parameters required to build a Kynema-FMB turbine
struct TurbineSetup
{
    std::string label;
    std::string input_file;
    std::string checkpoint_file;
    std::string controller_shared_lib_path;
    std::string controller_input_file;

    int num_blades;
    int num_blade_elem;
    int num_tower_elem;
    int num_pts_blade;
    int num_pts_tower;
    int section_refinement_blade;
    int section_refinement_tower;

    amrex::Real dt_ext;
    amrex::Real rotational_speed;
    amrex::Real generator_power;
    amrex::Real wind_speed;
    amrex::Real yaw;
    amrex::Real generator_efficiency;

    amrex::Real damping_factor;
    int nl_iter_max;
    amrex::Real abs_err_tol;
    amrex::Real rel_err_tol;
    amrex::Real gravity[3];
};

//! Implemented in KynemaFMBHost.cpp
TurbineInterfacePtr create_turbine_interface(const TurbineSetup& setup);

void set_time_step_index(TurbineInterface& ti, int time_index);

//! Number of actuator points (hub + aerodynamic nodes)
int num_force_points(TurbineInterface& ti);

void get_hub_state(
    TurbineInterface& ti, float* hub_abs_pos, float* hub_rot_vel);

void write_checkpoint_file(TurbineInterface& ti, const std::string& fname);

void calculate_motion(TurbineInterface& ti);

void calculate_aerodynamic_loads(TurbineInterface& ti, amrex::Real density);

void calculate_nodal_loads(TurbineInterface& ti);

bool step(TurbineInterface& ti);

void apply_controller(TurbineInterface& ti, amrex::Real time);

void write_output(TurbineInterface& ti);

/** Buffer exchange; each vector buffer stores the x, y, z components of all
 *  npts points one after the other, orientation stores 9 entries per point
 */
void populate_buffers(
    TurbineInterface& ti,
    int npts,
    bool one_way_coupled,
    float* chord,
    float* pos_force,
    float* pos_vel,
    float* solid_vel,
    float* force,
    float* orientation);

void pass_fluid_velocity_and_hub_load(
    TurbineInterface& ti, int npts, const float* fluid_vel, const float* force);

void save_and_zero_aerodynamic_loads(
    TurbineInterface& ti, int npts, float* force);

//! Implemented in KynemaFMBIface.cpp so that the host code can use AMReX
[[noreturn]] void abort_run(const std::string& msg);

void print_message(const std::string& msg);

amrex::Real interp_linear(
    const amrex::Real* xinp,
    const amrex::Real* yinp,
    int n,
    amrex::Real xout);

} // namespace sgf_fmb

#endif /* KYNEMA_FMB_HOST_H */
