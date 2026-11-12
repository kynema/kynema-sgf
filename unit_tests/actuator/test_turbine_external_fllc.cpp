#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <vector>
#include "ks_test_utils/MeshTest.H"
#include "AMReX_ParallelDescriptor.H"

// The external turbine data hold an MPI communicator; without MPI, AMReX
// provides stand-ins
#ifndef BL_USE_MPI
using amrex::MPI_Comm;
using amrex::MPI_COMM_NULL;
#endif

#include "src/actuator/turbine/external/ExtTurbIface.H"
#include "src/actuator/turbine/external/TurbineExternal.H"
#include "src/actuator/turbine/external/turbine_external_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

/** Turbine of a stand-in external solver whose blades bend
 *
 *  The solver gives the two point sets of a turbine as OpenFAST does through
 *  ExternalInflow, each list starting with the hub, then the blades from root
 *  to tip: the force points and the velocity points (the AeroDyn nodes). At
 *  every step the solver puts both sets of a blade on the same line, the one
 *  of that step. The blades bend out of the rotor plane a little more at
 *  every step and keep their length, so the straight-line distance of a point
 *  from the root of its blade changes from one step to the next.
 */
struct BendingTurbine : public ::ext_turb::ExternalTurbine
{
    float* position_at_vel(int dir) const override
    {
        return vel_pos[dir].data();
    }
    float* solid_velocity(int dir) const override
    {
        return solid_vel[dir].data();
    }
    float* fluid_velocity(int dir) const override
    {
        return fluid_vel[dir].data();
    }
    float* force(int dir) const override { return point_force[dir].data(); }
    float* position_at_force(int dir) const override
    {
        return force_pos[dir].data();
    }
    float* chord_at_force() const override { return chord.data(); }
    float* orientation() const override { return orient.data(); }

    int length_fluid_velocity(int /*dir*/) const override
    {
        return num_vel_pts();
    }
    int length_force(int /*dir*/) const override { return num_force_pts(); }
    int length_position_at_force(int /*dir*/) const override
    {
        return num_force_pts();
    }
    int length_orientation() const override { return 9 * num_force_pts(); }
    int num_vel_pts_blade() const override { return model_num_vel_pts_blade; }

    //! The hub and the force points of the blades
    int num_force_pts() const
    {
        return 1 + model_num_blades * model_num_pts_blade;
    }

    //! The hub and the velocity points of the blades
    int num_vel_pts() const
    {
        return 1 + model_num_blades * model_num_vel_pts_blade;
    }

    //! Curvature of the blades at a step of the solver
    amrex::Real curvature(const int step) const
    {
        return static_cast<amrex::Real>(step) * curvature_per_step;
    }

    //! Arc length from the root of point k of npts points evenly spaced
    //! along a blade
    amrex::Real arc_length(const int k, const int npts) const
    {
        return blade_length * static_cast<amrex::Real>(k) /
               static_cast<amrex::Real>(npts - 1);
    }

    /** Position of a point of a blade at a step of the solver, relative to
     *  the base of the turbine
     *
     *  Blade ib leaves the hub in the rotor plane (x = 0) at an azimuth of
     *  2 pi ib / num_blades. It is bent downwind into a circular arc that
     *  starts at its root, so a point at arc length s from the root is at a
     *  straight-line distance of 2 sin(kappa s / 2) / kappa from it.
     */
    std::array<float, 3>
    blade_point(const int ib, const amrex::Real s, const int step) const
    {
        const amrex::Real kappa = curvature(step);
        const amrex::Real azimuth = 2.0_rt * std::numbers::pi_v<amrex::Real> *
                                    static_cast<amrex::Real>(ib) /
                                    static_cast<amrex::Real>(model_num_blades);
        const amrex::Real radial =
            root_radius + ((kappa > 0.0_rt) ? std::sin(kappa * s) / kappa : s);
        const amrex::Real downwind =
            (kappa > 0.0_rt) ? (1.0_rt - std::cos(kappa * s)) / kappa : 0.0_rt;
        return {
            static_cast<float>(downwind),
            static_cast<float>(radial * std::sin(azimuth)),
            static_cast<float>(hub_height + radial * std::cos(azimuth))};
    }

    //! Buffers of a new turbine. Every point is oriented as the global frame.
    void allocate()
    {
        for (int dir = 0; dir < 3; ++dir) {
            force_pos[dir].assign(num_force_pts(), 0.0F);
            point_force[dir].assign(num_force_pts(), 0.0F);
            solid_vel[dir].assign(num_force_pts(), 0.0F);
            vel_pos[dir].assign(num_vel_pts(), 0.0F);
            fluid_vel[dir].assign(num_vel_pts(), 0.0F);
        }
        chord.assign(num_force_pts(), 1.0F);
        const auto npts = static_cast<std::size_t>(num_force_pts());
        orient.assign(9 * npts, 0.0F);
        for (std::size_t i = 0; i < npts; ++i) {
            orient[9 * i] = 1.0F;
            orient[(9 * i) + 4] = 1.0F;
            orient[(9 * i) + 8] = 1.0F;
        }
    }

    //! Both point sets at a step of the solver. With as many velocity points
    //! as force points on a blade, the two sets are the same points.
    void place_points(const int step)
    {
        const std::array<float, 3> hub{
            0.0F, 0.0F, static_cast<float>(hub_height)};
        set_point(force_pos, 0, hub);
        set_point(vel_pos, 0, hub);
        for (int ib = 0; ib < model_num_blades; ++ib) {
            for (int k = 0; k < model_num_pts_blade; ++k) {
                set_point(
                    force_pos, 1 + (ib * model_num_pts_blade) + k,
                    blade_point(ib, arc_length(k, model_num_pts_blade), step));
            }
            for (int k = 0; k < model_num_vel_pts_blade; ++k) {
                set_point(
                    vel_pos, 1 + (ib * model_num_vel_pts_blade) + k,
                    blade_point(
                        ib, arc_length(k, model_num_vel_pts_blade), step));
            }
        }
    }

    static void set_point(
        std::array<std::vector<float>, 3>& points,
        const int i,
        const std::array<float, 3>& pt)
    {
        for (int dir = 0; dir < 3; ++dir) {
            points[dir][i] = pt[dir];
        }
    }

    //! Blades of the model of the solver and points on each
    int model_num_blades{3};
    int model_num_pts_blade{5};
    int model_num_vel_pts_blade{5};

    //! Height of the hub above the base, distance of the roots of the blades
    //! from the hub, length of the blades
    amrex::Real hub_height{10.0_rt};
    amrex::Real root_radius{1.0_rt};
    amrex::Real blade_length{4.0_rt};

    //! Curvature the blades gain at every step of the solver
    amrex::Real curvature_per_step{0.05_rt};

    //! Buffers exchanged with the CFD solver, one list per direction
    mutable std::array<std::vector<float>, 3> force_pos;
    mutable std::array<std::vector<float>, 3> vel_pos;
    mutable std::array<std::vector<float>, 3> point_force;
    mutable std::array<std::vector<float>, 3> fluid_vel;
    mutable std::array<std::vector<float>, 3> solid_vel;
    mutable std::vector<float> chord;
    mutable std::vector<float> orient;
};

struct BendingSolverData
{};

} // namespace kynema_sgf_tests

namespace ext_turb {

template <>
std::string ext_id<::kynema_sgf_tests::BendingSolverData>()
{
    return "BendingSolver";
}

template <>
ExtTurbIface<
    ::kynema_sgf_tests::BendingTurbine,
    ::kynema_sgf_tests::BendingSolverData>::~ExtTurbIface() = default;

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::BendingTurbine,
    ::kynema_sgf_tests::BendingSolverData>::
    parse_inputs(
        const kynema_sgf::CFDSim& /*sim*/, const std::string& /*inp_name*/)
{}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::BendingTurbine,
    ::kynema_sgf_tests::BendingSolverData>::allocate_ext_turbines()
{
    m_is_initialized = true;
}

// The initial solution leaves the blades as they are at step 0, straight
template <>
void ExtTurbIface<
    ::kynema_sgf_tests::BendingTurbine,
    ::kynema_sgf_tests::BendingSolverData>::init_solution(const int local_id)
{
    auto& fi = *m_turbine_data[local_id];
    fi.place_points(fi.time_index);
    fi.is_solution0 = false;
}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::BendingTurbine,
    ::kynema_sgf_tests::BendingSolverData>::
    get_hub_stats(const int /*local_id*/)
{}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::BendingTurbine,
    ::kynema_sgf_tests::BendingSolverData>::
    calc_aero_loads(::kynema_sgf_tests::BendingTurbine& /*fi*/)
{}

// A step bends the blades further. As in OpenFAST, the step leaves both point
// sets at the new time: ExternalInflow sets the positions of the velocity
// points and of the force points together (SetExtInfwPositions).
template <>
void ExtTurbIface<
    ::kynema_sgf_tests::BendingTurbine,
    ::kynema_sgf_tests::BendingSolverData>::
    do_turbine_step(::kynema_sgf_tests::BendingTurbine& fi)
{
    fi.place_points(fi.time_index + 1);
}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::BendingTurbine,
    ::kynema_sgf_tests::BendingSolverData>::
    write_turbine_checkpoint(int& /*tid*/)
{}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::BendingTurbine,
    ::kynema_sgf_tests::BendingSolverData>::
    ext_init_turbine(::kynema_sgf_tests::BendingTurbine& fi)
{
    fi.time_index = 0;
    // The model of the solver gives the number of blades, as OpenFAST does
    fi.num_blades = fi.model_num_blades;
    fi.allocate();
    fi.place_points(0);
}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::BendingTurbine,
    ::kynema_sgf_tests::BendingSolverData>::
    ext_restart_turbine(::kynema_sgf_tests::BendingTurbine& /*fi*/)
{}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::BendingTurbine,
    ::kynema_sgf_tests::BendingSolverData>::
    ext_replay_turbine(::kynema_sgf_tests::BendingTurbine& /*fi*/)
{}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::BendingTurbine,
    ::kynema_sgf_tests::BendingSolverData>::
    prepare_netcdf_file(::kynema_sgf_tests::BendingTurbine& /*fi*/)
{}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::BendingTurbine,
    ::kynema_sgf_tests::BendingSolverData>::
    write_velocity_data(const ::kynema_sgf_tests::BendingTurbine& /*fi*/)
{}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::BendingTurbine,
    ::kynema_sgf_tests::BendingSolverData>::
    read_velocity_data(
        ::kynema_sgf_tests::BendingTurbine& /*fi*/,
        const ncutils::NCFile& /*ncf*/,
        const size_t /*tid*/)
{}

template class ExtTurbIface<
    ::kynema_sgf_tests::BendingTurbine,
    ::kynema_sgf_tests::BendingSolverData>;

} // namespace ext_turb

namespace kynema_sgf_tests {
namespace {

using BendingIface =
    ::ext_turb::ExtTurbIface<BendingTurbine, BendingSolverData>;
// The data of an actuator of this solver, as TurbineFast is for OpenFAST
using BendingData = kynema_sgf::actuator::
    TurbineFromExtSolver<BendingTurbine, BendingSolverData>::DataType;

//! Straight-line distance from the root of a point at arc length s along a
//! blade bent into a circular arc of curvature kappa
amrex::Real distance_from_root(const amrex::Real kappa, const amrex::Real s)
{
    return (kappa > 0.0_rt) ? 2.0_rt * std::sin(0.5_rt * kappa * s) / kappa : s;
}

class TurbineExternalFLLCTest : public MeshTest
{
protected:
    /** Turbine of three blades with the filtered lifting line correction, on
     *  the root process, initialised as Actuator::post_init_actions does it
     *
     *  \param num_pts_blade Force points on each blade
     *  \param num_vel_pts_blade Velocity points on each blade
     */
    void setup_turbine(
        BendingData& data,
        BendingIface& iface,
        const int num_pts_blade,
        const int num_vel_pts_blade)
    {
        namespace vs = kynema_sgf::vs;
        auto& info = data.info();
        auto& meta = data.meta();
        auto& fi = meta.ext_data;
        fi.model_num_blades = m_num_blades;
        fi.model_num_pts_blade = num_pts_blade;
        fi.model_num_vel_pts_blade = num_vel_pts_blade;

        const int iproc = amrex::ParallelDescriptor::MyProc();
        info.root_proc = iproc;
        info.is_root_proc = true;
        info.actuator_in_proc = true;
        info.procs = {iproc};

        // Base of the turbine, and the box that the inputs give it
        // (utils::read_inputs): 1.25 rotor radii around the base, up to the
        // hub height above it. The distances along a blade do not depend on
        // the base. At the origin the positions lose little to rounding in a
        // build in single precision.
        const vs::Vector base = vs::Vector::zero();
        const amrex::Real reach = 1.25_rt * (fi.root_radius + fi.blade_length);
        info.base_pos = base;
        info.bound_box = amrex::RealBox(
            base.x() - reach, base.y() - reach, base.z() - reach,
            base.x() + reach, base.y() + reach,
            base.z() + reach + fi.hub_height);

        // What the inputs leave: points of each blade, epsilon_chord, and
        // the data of the correction of every blade. The correction starts
        // at the first force call past fllc_start_time.
        meta.num_blades = m_num_blades;
        meta.num_pts_blade = num_pts_blade;
        meta.num_pts_tower = 0;
        meta.eps_chord = vs::Vector{1.0_rt, 1.0_rt, 1.0_rt};
        kynema_sgf::actuator::FLLCData fllc_inputs;
        fllc_inputs.fllc_start_time = 0.05_rt;
        meta.fllc.assign(m_num_blades, fllc_inputs);

        fi.tlabel = "T0";
        fi.tid_global = 0;
        fi.sim_mode = ::ext_turb::SimMode::init;
        fi.num_pts_tower = 0;
        fi.dt_cfd = 0.1_rt;
        fi.dt_ext = 0.1_rt;
        fi.num_substeps = 1;
        meta.ext_ptr = &iface;
        iface.register_turbine(fi);

        kynema_sgf::actuator::ops::init_data_op(data);
    }

    //! What Actuator::post_init_actions does for the turbine: take the
    //! velocity points from the solver, hand it the velocity sampled there,
    //! let it compute its loads and take its forces
    static void post_init(BendingData& data)
    {
        kynema_sgf::actuator::ops::update_pos_op(data);
        kynema_sgf::actuator::ops::update_fields_op(data);
        kynema_sgf::actuator::external::ext_loads(data);
        communicate(data);
    }

    //! What Actuator::pre_advance_work does for the turbine: as
    //! post_init_actions, with the step of the solver before its forces are
    //! taken
    static void pre_advance(BendingData& data)
    {
        kynema_sgf::actuator::ops::update_pos_op(data);
        kynema_sgf::actuator::ops::update_fields_op(data);
        kynema_sgf::actuator::external::ext_loads(data);
        kynema_sgf::actuator::external::ext_step(data);
        communicate(data);
    }

    static void communicate(BendingData& data)
    {
        kynema_sgf::actuator::ops::CommunicateModelDataOp<
            kynema_sgf::actuator::TurbineFromExtSolver<
                BendingTurbine, BendingSolverData>,
            kynema_sgf::actuator::ActSrcLine>()(data);
    }

    //! A CFD time step, as far as the turbine goes
    void time_step(BendingData& data)
    {
        auto& time = sim().time();
        time.new_timestep();
        time.set_current_cfl(0.25_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        pre_advance(data);
    }

    /** Run the turbine up to the start of the correction
     *
     *  The force calls at the initial time and in the first time step (at
     *  time 0) come before fllc_start_time. The correction starts in the
     *  second time step (at time 0.1), whose velocity points were taken from
     *  the solver at its step 1, before it took step 2.
     */
    void run_to_start_of_correction(BendingData& data)
    {
        post_init(data);
        time_step(data);
        for (const auto& fllc : data.meta().fllc) {
            ASSERT_FALSE(fllc.initialized);
        }
        time_step(data);
        ASSERT_EQ(data.meta().ext_data.time_index, 2);
        for (const auto& fllc : data.meta().fllc) {
            ASSERT_TRUE(fllc.initialized);
        }
    }

    const int m_num_blades{3};

    //! The solver gives the positions in single precision: the distances
    //! along a blade are off by up to about 5e-7
    const amrex::Real m_tol{1.0e-5_rt};
};

} // namespace

// fllc_init tells whether the velocity points of a blade are its force
// points, from the distance of every point from the first one of its set,
// and keeps these distances for the rest of the run. A blade that bends has
// other distances a step later, so both sets have to be those of one step:
// the force points of the step the solver has just taken, and its velocity
// points of that step.
TEST_F(TurbineExternalFLLCTest, same_points_of_bending_blade_stay_same_points)
{
    const int npb = 5;
    initialize_mesh();
    BendingIface iface(sim());
    BendingData data(sim(), "T0", 0);
    setup_turbine(data, iface, npb, npb);
    const auto& fi = data.meta().ext_data;

    // From step 1 to step 2 the tip of a blade comes 0.02 closer to its
    // root, 25 times the tolerance with which fllc_init tells two sets of
    // distances apart: 0.001 of the distance of the last point over the
    // number of points
    const amrex::Real tip = fi.arc_length(npb - 1, npb);
    const amrex::Real tip_distance = distance_from_root(fi.curvature(2), tip);
    ASSERT_GT(
        distance_from_root(fi.curvature(1), tip) - tip_distance,
        10.0_rt * 1.0e-3_rt * tip_distance / static_cast<amrex::Real>(npb));

    ASSERT_NO_FATAL_FAILURE(run_to_start_of_correction(data));

    // Five AeroDyn nodes and five force points at the nodes, on every blade:
    // the same points, which the correction is handed to one by one
    const auto kappa = fi.curvature(2);
    for (int ib = 0; ib < m_num_blades; ++ib) {
        const auto& fllc = data.meta().fllc[ib];
        EXPECT_FALSE(fllc.mismatched_points) << "blade " << ib;
        ASSERT_EQ(static_cast<int>(fllc.span_distance_force.size()), npb);
        ASSERT_EQ(static_cast<int>(fllc.span_distance_vel.size()), npb);
        for (int k = 0; k < npb; ++k) {
            const auto gold = distance_from_root(kappa, fi.arc_length(k, npb));
            EXPECT_NEAR(fllc.span_distance_force[k], gold, m_tol)
                << "blade " << ib << ", force point " << k;
            EXPECT_NEAR(fllc.span_distance_vel[k], gold, m_tol)
                << "blade " << ib << ", velocity point " << k;
        }
    }
}

// Seven velocity points and five force points along the bending blade. The
// correction interpolates between the two sets by the distances fllc_init
// keeps, for the rest of the run, so both sets of distances have to be those
// of the step of the force points.
TEST_F(TurbineExternalFLLCTest, distances_of_other_points_are_of_same_step)
{
    const int npb = 5;
    const int nvb = 7;
    initialize_mesh();
    BendingIface iface(sim());
    BendingData data(sim(), "T0", 0);
    setup_turbine(data, iface, npb, nvb);
    const auto& fi = data.meta().ext_data;

    ASSERT_NO_FATAL_FAILURE(run_to_start_of_correction(data));

    const auto kappa = fi.curvature(2);
    for (int ib = 0; ib < m_num_blades; ++ib) {
        const auto& fllc = data.meta().fllc[ib];
        EXPECT_TRUE(fllc.mismatched_points) << "blade " << ib;
        ASSERT_EQ(static_cast<int>(fllc.span_distance_force.size()), npb);
        ASSERT_EQ(static_cast<int>(fllc.span_distance_vel.size()), nvb);
        for (int k = 0; k < npb; ++k) {
            EXPECT_NEAR(
                fllc.span_distance_force[k],
                distance_from_root(kappa, fi.arc_length(k, npb)), m_tol)
                << "blade " << ib << ", force point " << k;
        }
        for (int j = 0; j < nvb; ++j) {
            EXPECT_NEAR(
                fllc.span_distance_vel[j],
                distance_from_root(kappa, fi.arc_length(j, nvb)), m_tol)
                << "blade " << ib << ", velocity point " << j;
        }
    }
}

} // namespace kynema_sgf_tests
