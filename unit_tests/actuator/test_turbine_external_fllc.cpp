#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <numbers>
#include <string>
#include <vector>
#include "ks_test_utils/MeshTest.H"
#include "AMReX_ParallelDescriptor.H"

// MPI communicator in the turbine data; AMReX stand-ins without MPI
#ifndef BL_USE_MPI
using amrex::MPI_Comm;
using amrex::MPI_COMM_NULL;
#endif

#include "src/actuator/turbine/external/ExtTurbIface.H"
#include "src/actuator/turbine/external/TurbineExternal.H"
#include "src/actuator/turbine/external/turbine_external_ops.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"
#include "AMReX_Utility.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

/** Stand-in external solver with bending blades
 *
 *  Force and velocity points (AeroDyn nodes) as OpenFAST gives them through
 *  ExternalInflow: hub first, then blades root to tip; both sets on the line
 *  of the current step. Blades bend further out of plane each step at fixed
 *  length, so straight-line distances from the root change per step.
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
        return 1 + (model_num_blades * model_num_pts_blade);
    }

    //! The hub and the velocity points of the blades
    int num_vel_pts() const
    {
        return 1 + (model_num_blades * model_num_vel_pts_blade);
    }

    //! Curvature of the blades at a step of the solver
    amrex::Real curvature(const int step) const
    {
        return static_cast<amrex::Real>(step) * curvature_per_step;
    }

    //! Arc length from the root of point k of npts evenly spaced points
    amrex::Real arc_length(const int k, const int npts) const
    {
        return blade_length * static_cast<amrex::Real>(k) /
               static_cast<amrex::Real>(npts - 1);
    }

    /** Blade point position at a solver step, relative to the turbine base
     *
     *  Blade ib: in the rotor plane (x = 0) at azimuth 2 pi ib / num_blades,
     *  bent downwind into a circular arc from its root; arc length s is at
     *  straight-line distance 2 sin(kappa s / 2) / kappa from the root.
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
            static_cast<float>(hub_height + (radial * std::cos(azimuth)))};
    }

    //! Buffers of a new turbine; points oriented as the global frame
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

    //! Both point sets at a solver step (identical with equal counts)
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

    //! Hub height above the base, root distance from the hub, blade length
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

// A step bends the blades further and, as in OpenFAST (SetExtInfwPositions),
// moves both point sets to the new time.
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
    //! New interface, heap-owned by the fixture (view data point to it)
    BendingIface& make_iface()
    {
        m_ifaces.push_back(std::make_unique<BendingIface>(sim()));
        return *m_ifaces.back();
    }

    /** Three-blade turbine with the FLLC on the root process, initialized as
     *  by Actuator::post_init_actions
     *
     *  \param num_pts_blade Force points on each blade
     *  \param num_vel_pts_blade Velocity points on each blade
     *  \param sim_mode How the solver starts the turbine
     */
    void setup_turbine(
        BendingData& data,
        BendingIface& iface,
        const int num_pts_blade,
        const int num_vel_pts_blade,
        const ::ext_turb::SimMode sim_mode = ::ext_turb::SimMode::init) const
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

        // Base and its box (utils::read_inputs): 1.25 rotor radii around the
        // base, up to hub height. Base at the origin to limit single-precision
        // rounding; distances along a blade do not depend on it.
        const vs::Vector base = vs::Vector::zero();
        const amrex::Real reach = 1.25_rt * (fi.root_radius + fi.blade_length);
        info.base_pos = base;
        info.bound_box = amrex::RealBox(
            base.x() - reach, base.y() - reach, base.z() - reach,
            base.x() + reach, base.y() + reach,
            base.z() + reach + fi.hub_height);

        // From the inputs: points per blade, epsilon_chord, per-blade FLLC
        // data. FLLC starts at the first force call past fllc_start_time.
        meta.num_blades = m_num_blades;
        meta.num_pts_blade = num_pts_blade;
        meta.num_pts_tower = 0;
        meta.eps_chord = vs::Vector{1.0_rt, 1.0_rt, 1.0_rt};
        kynema_sgf::actuator::FLLCData fllc_inputs;
        fllc_inputs.fllc_start_time = 0.05_rt;
        meta.fllc.assign(m_num_blades, fllc_inputs);

        fi.tlabel = "T0";
        fi.tid_global = 0;
        fi.sim_mode = sim_mode;
        fi.num_pts_tower = 0;
        fi.dt_cfd = 0.1_rt;
        fi.dt_ext = 0.1_rt;
        fi.num_substeps = 1;
        meta.ext_ptr = &iface;
        iface.register_turbine(fi);

        kynema_sgf::actuator::ops::init_data_op(data);
    }

    //! As Actuator::post_init_actions: velocity points from the solver,
    //! sampled velocity to it, loads, forces back
    static void post_init(BendingData& data)
    {
        kynema_sgf::actuator::ops::update_pos_op(data);
        kynema_sgf::actuator::ops::update_fields_op(data);
        kynema_sgf::actuator::external::ext_loads(data);
        communicate(data);
    }

    //! As Actuator::pre_advance_work: post_init_actions plus a solver step
    //! before the forces
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

    /** Run the turbine up to the FLLC start
     *
     *  Force calls at init and in step 1 (time 0) precede fllc_start_time; the
     *  FLLC starts in step 2 (time 0.1), with velocity points from solver step
     *  1, taken before solver step 2.
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

    /** Run to the FLLC start, mark each blade's FLLC state and forces
     *  (m_force_mark) and write them into checkpoint chkname, as the
     *  end-of-step CFD checkpoint does
     *
     *  \param chk_time Time of the checkpoint
     */
    void write_marked_state(
        const std::string& chkname,
        const int num_pts_blade,
        const kynema_sgf::vs::Vector& mark,
        amrex::Real& chk_time)
    {
        auto& iface = make_iface();
        BendingData data(sim(), "T0", 0);
        setup_turbine(data, iface, num_pts_blade, num_pts_blade);
        ASSERT_NO_FATAL_FAILURE(run_to_start_of_correction(data));
        for (auto& fllc : data.meta().fllc) {
            fllc.correction_velocity.assign(
                fllc.correction_velocity.size(), mark);
        }
        for (auto& blade : data.meta().blades) {
            for (auto& f : blade.force) {
                f = m_force_mark;
            }
        }
        amrex::UtilCreateDirectory(chkname, 0755);
        kynema_sgf::actuator::ops::external_write_checkpoint(data, chkname);
        chk_time = sim().time().new_time();
        ASSERT_TRUE(
            kynema_sgf::actuator::restart::has_file(chkname, "T0", chk_time));
    }

    //! Restart settings from checkpoint chkname before the first force call:
    //! checkpoint name and time
    void restart_flow_from(const std::string& chkname, const amrex::Real time)
    {
        amrex::ParmParse pp("io");
        pp.add("restart_file", chkname);
        sim().io_manager().initialize_io();
        sim().time().set_restart_time(sim().time().time_index(), time);
    }

    const int m_num_blades{3};

    //! Blade force written with the marked state
    const kynema_sgf::vs::Vector m_force_mark{-1.5_rt, 0.5_rt, -0.25_rt};

    std::vector<std::unique_ptr<BendingIface>> m_ifaces;

    //! Single-precision solver positions: distances off by up to ~5e-7
    const amrex::Real m_tol{1.0e-5_rt};
};

} // namespace

// fllc_init compares and keeps point distances from each set's first
// point. A bending blade changes them per step, so both sets must be of the
// step the solver just took.
TEST_F(TurbineExternalFLLCTest, same_points_of_bending_blade_stay_same_points)
{
    // Serial: the turbine data would make every rank its own root
    if (amrex::ParallelDescriptor::NProcs() > 1) {
        GTEST_SKIP();
    }
    const int npb = 5;
    initialize_mesh();
    auto& iface = make_iface();
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

    // Five AeroDyn nodes and five force points at them, per blade: one set
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

// Seven velocity and five force points on a bending blade. The FLLC
// interpolates by the distances fllc_init keeps, so both must be of the
// force points' step.
TEST_F(TurbineExternalFLLCTest, distances_of_other_points_are_of_same_step)
{
    // Serial: the turbine data would make every rank its own root
    if (amrex::ParallelDescriptor::NProcs() > 1) {
        GTEST_SKIP();
    }
    const int npb = 5;
    const int nvb = 7;
    initialize_mesh();
    auto& iface = make_iface();
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

// Restart from a checkpoint with FLLC state, solver at the checkpoint time
// (replay): the FLLC continues from that state at the first force call.
TEST_F(TurbineExternalFLLCTest, restart_goes_on_from_state_in_checkpoint)
{
    // Serial: the turbine data would make every rank its own root
    if (amrex::ParallelDescriptor::NProcs() > 1) {
        GTEST_SKIP();
    }
    namespace vs = kynema_sgf::vs;
    const int npb = 5;
    const std::string chkname = "test_turbine_external_fllc_replay_chk";
    const vs::Vector mark{0.25_rt, -0.5_rt, 0.75_rt};
    initialize_mesh();
    amrex::Real chk_time = 0.0_rt;
    ASSERT_NO_FATAL_FAILURE(write_marked_state(chkname, npb, mark, chk_time));
    restart_flow_from(chkname, chk_time);

    auto& iface = make_iface();
    BendingData data(sim(), "T0", 0);
    setup_turbine(data, iface, npb, npb, ::ext_turb::SimMode::replay);
    post_init(data);

    for (int ib = 0; ib < m_num_blades; ++ib) {
        const auto& fllc = data.meta().fllc[ib];
        ASSERT_TRUE(fllc.initialized) << "blade " << ib;
        for (const auto& u : fllc.correction_velocity) {
            EXPECT_EQ(u.x(), mark.x()) << "blade " << ib;
            EXPECT_EQ(u.y(), mark.y()) << "blade " << ib;
            EXPECT_EQ(u.z(), mark.z()) << "blade " << ib;
        }
        // Forces of the checkpoint, not those the solver gives at restart
        for (const auto& f : data.meta().blades[ib].force) {
            EXPECT_EQ(f.x(), m_force_mark.x()) << "blade " << ib;
            EXPECT_EQ(f.y(), m_force_mark.y()) << "blade " << ib;
            EXPECT_EQ(f.z(), m_force_mark.z()) << "blade " << ib;
        }
    }
}

// Same restart, solver starting the turbine anew (init): state unused; the
// FLLC starts from zero in step 1, as without a state.
TEST_F(TurbineExternalFLLCTest, restart_of_turbine_started_anew_ignores_state)
{
    // Serial: the turbine data would make every rank its own root
    if (amrex::ParallelDescriptor::NProcs() > 1) {
        GTEST_SKIP();
    }
    namespace vs = kynema_sgf::vs;
    const int npb = 5;
    const std::string chkname = "test_turbine_external_fllc_init_chk";
    const vs::Vector mark{0.25_rt, -0.5_rt, 0.75_rt};
    initialize_mesh();
    amrex::Real chk_time = 0.0_rt;
    ASSERT_NO_FATAL_FAILURE(write_marked_state(chkname, npb, mark, chk_time));
    restart_flow_from(chkname, chk_time);

    auto& iface = make_iface();
    BendingData data(sim(), "T0", 0);
    setup_turbine(data, iface, npb, npb, ::ext_turb::SimMode::init);
    post_init(data);
    for (int ib = 0; ib < m_num_blades; ++ib) {
        EXPECT_FALSE(data.meta().fllc[ib].initialized) << "blade " << ib;
    }

    time_step(data);
    for (int ib = 0; ib < m_num_blades; ++ib) {
        const auto& fllc = data.meta().fllc[ib];
        ASSERT_TRUE(fllc.initialized) << "blade " << ib;
        for (const auto& u : fllc.correction_velocity) {
            EXPECT_EQ(u.x(), 0.0_rt) << "blade " << ib;
            EXPECT_EQ(u.y(), 0.0_rt) << "blade " << ib;
            EXPECT_EQ(u.z(), 0.0_rt) << "blade " << ib;
        }
    }
}

} // namespace kynema_sgf_tests
