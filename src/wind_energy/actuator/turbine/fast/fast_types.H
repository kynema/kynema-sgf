#ifndef FAST_TYPES_H
#define FAST_TYPES_H

#include "amr-wind/wind_energy/actuator/turbine/fast/fast_wrapper.H"
#include "amr-wind/wind_energy/actuator/turbine/external/external_base_types.H"

#include <string>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace ext_turb {

/** Representation of a turbine for exchanging data with FAST
 */
struct FastTurbine : public ext_turb::ExternalTurbine
{
    /** Variables exclusive to OpenFAST turbine type
     */

    //! Node cluster type for the chord
    int chord_cluster_type{0};

    exw_fast::OfInputType to_cfd;
    exw_fast::OfOutputType from_cfd;

    float* position_at_vel(int dir) const override
    {
        return (
            dir == 0 ? to_cfd.pxVel : (dir == 1 ? to_cfd.pyVel : to_cfd.pzVel));
    }

    float* solid_velocity(int dir) const override
    {
        return (
            dir == 0 ? to_cfd.xdotForce
                     : (dir == 1 ? to_cfd.ydotForce : to_cfd.zdotForce));
    }

    float* fluid_velocity(int dir) const override
    {
        return (dir == 0 ? from_cfd.u : (dir == 1 ? from_cfd.v : from_cfd.w));
    }

    float* force(int dir) const override
    {
        return (dir == 0 ? to_cfd.fx : (dir == 1 ? to_cfd.fy : to_cfd.fz));
    }

    float* position_at_force(int dir) const override
    {
        return (
            dir == 0 ? to_cfd.pxForce
                     : (dir == 1 ? to_cfd.pyForce : to_cfd.pzForce));
    }

    float* chord_at_force() const override { return to_cfd.forceNodesChord; }

    float* orientation() const override { return to_cfd.pOrientation; }

    int length_fluid_velocity(int dir) const override
    {
        return (
            dir == 0 ? from_cfd.u_Len
                     : (dir == 1 ? from_cfd.v_Len : from_cfd.w_Len));
    }

    int length_force(int dir) const override
    {
        return (
            dir == 0 ? to_cfd.fx_Len
                     : (dir == 1 ? to_cfd.fy_Len : to_cfd.fz_Len));
    }

    int length_position_at_force(int dir) const override
    {
        return (
            dir == 0 ? to_cfd.pxForce_Len
                     : (dir == 1 ? to_cfd.pyForce_Len : to_cfd.pzForce_Len));
    }

    int length_orientation() const override { return to_cfd.pOrientation_Len; }

    int num_vel_pts_blade() const override { return num_blade_elem; }

#if OPENFAST_VERSION_MAJOR != 4
    exw_fast::SC_DX_InputType to_sc;
    exw_fast::SC_DX_OutputType from_sc;
#endif
};

struct FastSolverData
{
    std::string output_dir{"fast_velocity_data"};

    int m_num_sc_inputs{0};
    int m_num_sc_outputs{0};

#if OPENFAST_VERSION_MAJOR == 4
    int m_inflow_type{2};
#endif

    int m_num_sc_inputs_glob{0};
    float m_init_sc_inputs_glob{0.0_rt};
    float m_init_sc_inputs_turbine{0.0_rt};
};

} // namespace ext_turb

#endif /* FAST_TYPES_H */
