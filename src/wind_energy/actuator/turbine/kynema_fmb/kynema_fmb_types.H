#ifndef KYNEMA_FMB_TYPES_H
#define KYNEMA_FMB_TYPES_H

#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include "src/wind_energy/actuator/turbine/external/external_base_types.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

/** Kynema-FMB types are only used in KynemaFMBHost.cpp, which is compiled with
 *  the same compiler as Kokkos (see CMakeLists.txt), since nvcc and the host
 *  compiler lay out Kokkos::View differently. This header is also compiled by
 *  that file, so it must not include AMReX headers that contain GPU code.
 */
namespace kynema_fmb::interfaces {
class TurbineInterface;
} // namespace kynema_fmb::interfaces

namespace sgf_fmb {

static std::string checkpoint_filename(int time_index, std::string turbine_name)
{
    std::stringstream idx_str;
    idx_str << "kynema_fmb_" << std::setfill('0') << std::setw(5) << time_index
            << "_" << turbine_name << ".chk";
    return idx_str.str();
}
} // namespace sgf_fmb

namespace ext_turb {

/** Representation of a turbine for exchanging data with Kynema-FMB
 */
struct KynemaFMBTurbine : public ext_turb::ExternalTurbine
{
    /** Variables exclusive to KynemaFMB turbine type
     */
    std::unique_ptr<kynema_fmb::interfaces::TurbineInterface> interface;

    KynemaFMBTurbine();
    ~KynemaFMBTurbine();

    // Need buffers to make data contiguous to fit into coupling
    void allocate_buffers()
    {
        buffer_cf.resize(length_force(0));
        buffer_pv.resize(
            length_fluid_velocity(0) + length_fluid_velocity(1) +
            length_fluid_velocity(2));
        buffer_sv.resize(
            length_fluid_velocity(0) + length_fluid_velocity(1) +
            length_fluid_velocity(2));
        buffer_fv.resize(
            length_fluid_velocity(0) + length_fluid_velocity(1) +
            length_fluid_velocity(2));

        buffer_f.resize(length_force(0) + length_force(1) + length_force(2));
        buffer_pf.resize(
            length_position_at_force(0) + length_position_at_force(1) +
            length_position_at_force(2));
        buffer_ori.resize(length_orientation());

        ptr_cf = &buffer_cf[0];
        ptr_pv = &buffer_pv[0];
        ptr_sv = &buffer_sv[0];
        ptr_fv = &buffer_fv[0];
        ptr_f = &buffer_f[0];
        ptr_pf = &buffer_pf[0];
        ptr_ori = &buffer_ori[0];
    }

    void populate_buffers(bool one_way_coupled);

    void pass_fluid_velocity_and_hub_load();

    void save_and_zero_aerodynamic_loads();

    /** Helper functions for navigating indices */
    int body_index(int id) const;

    int node_index(int id) const;

    /** Data access functions needed for every turbine type */

    float* position_at_vel(int dir) const override
    {
        return ptr_pv + dir * length_fluid_velocity(0);
    }

    float* solid_velocity(int dir) const override
    {
        return ptr_sv + dir * length_fluid_velocity(0);
    }

    float* fluid_velocity(int dir) const override
    {
        return ptr_fv + dir * length_fluid_velocity(0);
    }

    float* force(int dir) const override
    {
        return ptr_f + dir * length_force(0);
    }

    float* position_at_force(int dir) const override
    {
        return ptr_pf + dir * length_position_at_force(0);
    }

    float* chord_at_force() const override { return ptr_cf; }

    float* orientation() const override { return ptr_ori; }

    int length_fluid_velocity(int dir) const override
    {
        return length_force(dir);
    }

    int length_force(int /*dir*/) const override;

    int length_position_at_force(int dir) const override
    {
        return length_force(dir);
    }

    int length_orientation() const override
    {
        return length_force(0) * AMREX_SPACEDIM * AMREX_SPACEDIM;
    }

    int num_vel_pts_blade() const override { return num_pts_blade; }

    amrex::Real fluid_density{1.0_rt};
    amrex::Real rotational_speed{0.0_rt};
    amrex::Real generator_power{0.0_rt};
    amrex::Real wind_speed{0.0_rt};
    amrex::Real yaw{0.0_rt};
    amrex::Real generator_efficiency{1.0_rt};
    int section_refinement_blade{0};
    int section_refinement_tower{2};
    int substep_counter{0};
    bool one_way_coupled{false};

    std::string controller_shared_lib_path{""};
    std::string controller_input_file{""};

    std::vector<float> buffer_cf;
    std::vector<float> buffer_pv;
    std::vector<float> buffer_sv;
    std::vector<float> buffer_fv;
    std::vector<float> buffer_f;
    std::vector<float> buffer_pf;
    std::vector<float> buffer_ori;
    float* ptr_cf{nullptr};
    float* ptr_pv{nullptr};
    float* ptr_sv{nullptr};
    float* ptr_fv{nullptr};
    float* ptr_f{nullptr};
    float* ptr_pf{nullptr};
    float* ptr_ori{nullptr};
};

struct KynemaFMBSolverData
{
    std::string output_dir{"kynema_fmb_velocity_data"};
    int m_offset_velocity{1};

    amrex::Real damping_factor{0.0_rt};
    int nl_iter_max{12};
    amrex::Real abs_err_tol{1.0e-5_rt};
    amrex::Real rel_err_tol{1.0e-4_rt};
    std::vector<amrex::Real> gravity{0.0_rt, 0.0_rt, -9.81_rt};
};

} // namespace ext_turb

//! Functions that use Kynema-FMB types, defined in KynemaFMBHost.cpp
namespace sgf_fmb {
void init_turbine_interface(
    ext_turb::KynemaFMBTurbine& fi,
    const ext_turb::KynemaFMBSolverData& solver_data);

void update_turbine(ext_turb::KynemaFMBTurbine& fi, bool advance);

void get_hub_state(ext_turb::KynemaFMBTurbine& fi);

void set_time_step_index(ext_turb::KynemaFMBTurbine& fi);

void write_checkpoint_file(
    ext_turb::KynemaFMBTurbine& fi, const std::string& fname);
} // namespace sgf_fmb

#endif /* KYNEMA_FMB_TYPES_H */
