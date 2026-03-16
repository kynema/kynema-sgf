#ifndef KYNEMA_TYPES_H
#define KYNEMA_TYPES_H

#include <string>
#include "amr-wind/wind_energy/actuator/turbine/kynema/kynema_wrapper.H"
#include "amr-wind/wind_energy/actuator/turbine/external/external_base_types.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace ext_turb {

/** Representation of a turbine for exchanging data with Kynema
 */
struct KynemaTurbine : public ext_turb::ExternalTurbine
{
    /** Variables exclusive to Kynema turbine type
     */
    std::unique_ptr<kynema::interfaces::TurbineInterface> interface;

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

    void populate_buffers()
    {
        // Get hub information for nacelle force
        // Chord is not used, set it to zero
        buffer_cf[0] = 0.;
        // Force is calculated by AMR-Wind, passed to Kynema
        // Get position directly from node information
        buffer_pf[0] = interface->Turbine().hub_node.position[0];
        buffer_pf[length_force(0)] = interface->Turbine().hub_node.position[1];
        buffer_pf[2 * length_force(0)] =
            interface->Turbine().hub_node.position[2];
        buffer_pv[0] = buffer_pf[0];
        buffer_pv[length_force(0)] = buffer_pf[length_force(0)];
        buffer_pv[2 * length_force(0)] = buffer_pf[2 * length_force(0)];
        // Get solid velocity as well, not used at the moment though
        buffer_sv[0] = interface->Turbine().hub_node.velocity[0];
        buffer_sv[length_force(0)] = interface->Turbine().hub_node.velocity[1];
        buffer_sv[2 * length_force(0)] =
            interface->Turbine().hub_node.velocity[2];

        // Get blade information
        int idx = 0;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            // Skip hub node
            ++idx;
            // Fluid velocity (fv) is not populated here, it comes from the
            // CFD solver. Need to add step to pass it from buffer to Kynema
            for (int n = 1; n < length_force(0); ++n) {
                if (dir == 0) {
                    buffer_cf[n] = (float)interface->Aerodynamics()
                                       .bodies[body_index(n)]
                                       .chord[node_index(n)];
                }
                buffer_f[idx] = (float)interface->Aerodynamics()
                                    .bodies[body_index(n)]
                                    .loads[node_index(n)][dir];
                buffer_pf[idx] = (float)interface->Aerodynamics()
                                     .bodies[body_index(n)]
                                     .x_motion[node_index(n)][dir];
                buffer_pv[idx] = buffer_pf[idx];
                buffer_sv[idx] = (float)interface->Aerodynamics()
                                     .bodies[body_index(n)]
                                     .v_motion[node_index(n)][dir];
                ++idx;
            }
        }
        // Orientation stores a 3x3 tensor (9 elements) for every point
        std::array<float, 4> q;
        for (int n = 0; n < length_force(0); ++n) {
            // convert quaternions to rotation matrix
            for (int iq = 0; iq < 4; ++iq) {
                if (n == 0) {
                    q[iq] =
                        (float)interface->Turbine().hub_node.position[iq + 3];
                } else {
                    q[iq] = (float)interface->Aerodynamics()
                                .bodies[body_index(n)]
                                .qqr_motion_map[node_index(n)][iq];
                }
            }
            const int idx = n * 9;
            buffer_ori[idx + 0] =
                q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
            buffer_ori[idx + 1] = 2. * (q[1] * q[2] - q[0] * q[3]);
            buffer_ori[idx + 2] = 2. * (q[1] * q[3] + q[0] * q[2]);
            buffer_ori[idx + 3] = 2. * (q[1] * q[2] + q[0] * q[3]);
            buffer_ori[idx + 4] =
                (q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3]);
            buffer_ori[idx + 5] = 2. * (q[2] * q[3] - q[0] * q[1]);
            buffer_ori[idx + 6] = 2. * (q[1] * q[3] - q[0] * q[2]);
            buffer_ori[idx + 7] = 2. * (q[2] * q[3] + q[0] * q[1]);
            buffer_ori[idx + 8] =
                (q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3]);
        }
    }

    void pass_fluid_velocity_and_hub_load()
    {
        // Fluid velocity at hub does not need to be passed, loads computed
        // internally in AMR-Wind
        int idx = 0;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            // Skip hub node
            ++idx;
            for (int n = 1; n < length_fluid_velocity(0); ++n) {
                interface->Aerodynamics()
                    .bodies[body_index(n)]
                    .v_inflow[node_index(n)][dir] = buffer_fv[idx];
                ++idx;
            }
        }
        // Hub load (also called nacelle_force) does need to be passed
        interface->Turbine().hub_node.loads[0] = buffer_f[0];
        interface->Turbine().hub_node.loads[1] = buffer_f[length_force(0)];
        interface->Turbine().hub_node.loads[2] =
            buffer_f[length_force(0) + length_force(1)];
        // Fluid velocity at hub is passed separately for the sake of controller
        std::array<amrex::Real, 3> hub_inflow{
            buffer_fv[0], buffer_fv[length_force(0)],
            buffer_fv[length_force(0) + length_force(1)]};
        interface->SetHubInflow(hub_inflow);
    }

    /** Helper functions for navigating indices */
    int body_index(int id) const
    {
        // Considering hub, which is 1 node, subtract
        id -= 1;
        int b_id = 0;
        int id_diff = id;
        for (int b = 0; b < (int)interface->Aerodynamics().bodies.size(); ++b) {
            id_diff -= static_cast<int>(
                interface->Aerodynamics().bodies[b].loads.size());
            if (id_diff < 0) {
                b_id = b;
                break;
            }
        }
        return b_id;
    }

    int node_index(int id) const
    {
        // Considering hub, which is 1 node, subtract
        id -= 1;
        int n_id = 0;
        int id_diff = id;
        for (int b = 0; b < (int)interface->Aerodynamics().bodies.size(); ++b) {
            id_diff -= (int)interface->Aerodynamics().bodies[b].loads.size();
            if (id_diff < 0) {
                n_id = id_diff +
                       (int)interface->Aerodynamics().bodies[b].loads.size();
                break;
            }
        }
        return n_id;
    }

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

    int length_force(int dir) const override
    {
        unsigned int n_sect{1};
        for (unsigned int n = 0; n < interface->Aerodynamics().bodies.size();
             ++n) {
            n_sect += interface->Aerodynamics().bodies[n].loads.size();
        }
        return static_cast<int>(n_sect);
    }

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
    int substep_counter{0};

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

struct KynemaSolverData
{
    std::string output_dir{"kynema_velocity_data"};

    amrex::Real damping_factor{0.0_rt};
    int nl_iter_max{12};
    amrex::Real abs_err_tol{1.0e-5_rt};
    amrex::Real rel_err_tol{1.0e-4_rt};
    amrex::Vector<amrex::Real> gravity{0.0_rt, 0.0_rt, -9.81_rt};
};

} // namespace ext_turb

#endif /* KYNEMA_TYPES_H */
