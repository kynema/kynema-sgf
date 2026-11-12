#include <array>
#include <string>
#include <vector>
#include "ks_test_utils/MeshTest.H"
#include "test_ext_turbine_utils.H"
#include "src/actuator/turbine/external/ExtTurbIface.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

/** Turbine of a stand-in external solver that counts what it is asked to do
 *
 *  The solver has no physics. An initial solution, a step, a checkpoint and a
 *  velocity record only add to the counters, so a test can tell which of them
 *  the interface requested and in which state it left the turbine.
 */
struct CountingTurbine : public ::ext_turb::ExternalTurbine
{
    float* position_at_vel(int /*dir*/) const override { return nullptr; }
    float* solid_velocity(int /*dir*/) const override { return nullptr; }
    float* fluid_velocity(int dir) const override { return &fluid_vel[dir]; }
    float* force(int dir) const override { return &hub_force[dir]; }
    float* position_at_force(int /*dir*/) const override { return nullptr; }
    float* chord_at_force() const override { return nullptr; }
    float* orientation() const override { return nullptr; }

    int length_fluid_velocity(int /*dir*/) const override { return 1; }
    int length_force(int /*dir*/) const override { return 1; }
    int length_position_at_force(int /*dir*/) const override { return 1; }
    int length_orientation() const override { return 0; }
    int num_vel_pts_blade() const override { return 0; }

    //! Step count of the solver at the time the CFD solver restarts from
    int time_index_at_restart{0};

    int num_initial_solutions{0};
    int num_steps{0};
    int num_checkpoints{0};

    //! CFD time step of each velocity record written
    mutable std::vector<int> velocity_records;

    //! Velocity at the hub and force of the hub, the only point of the turbine
    mutable std::array<float, 3> fluid_vel{0.0F, 0.0F, 0.0F};
    mutable std::array<float, 3> hub_force{0.0F, 0.0F, 0.0F};
};

struct CountingSolverData
{};

} // namespace kynema_sgf_tests

namespace ext_turb {

template <>
std::string ext_id<::kynema_sgf_tests::CountingSolverData>()
{
    return "CountingSolver";
}

template <>
ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::~ExtTurbIface() = default;

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::parse_inputs(
    const kynema_sgf::CFDSim& /*sim*/, const std::string& /*inp_name*/)
{}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::allocate_ext_turbines()
{
    m_is_initialized = true;
}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::init_solution(const int local_id)
{
    auto& fi = *m_turbine_data[local_id];
    ++fi.num_initial_solutions;
    fi.is_solution0 = false;
}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::get_hub_stats(
    const int /*local_id*/)
{}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::calc_aero_loads(
    ::kynema_sgf_tests::CountingTurbine& /*fi*/)
{}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::do_turbine_step(
    ::kynema_sgf_tests::CountingTurbine& fi)
{
    ++fi.num_steps;
}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::write_turbine_checkpoint(int& tid)
{
    ++m_turbine_data[tid]->num_checkpoints;
}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::ext_init_turbine(
    ::kynema_sgf_tests::CountingTurbine& fi)
{
    fi.time_index = 0;
}

// What reading a checkpoint of the solver does: the solver is at the step
// the checkpoint was written at
template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::ext_restart_turbine(
    ::kynema_sgf_tests::CountingTurbine& fi)
{
    fi.time_index = fi.time_index_at_restart;
}

// What a replay does: the initial solution, then the steps up to the restart
// time with the velocities of the earlier run
template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::ext_replay_turbine(
    ::kynema_sgf_tests::CountingTurbine& fi)
{
    ++fi.num_initial_solutions;
    fi.is_solution0 = false;
    for (; fi.time_index < fi.time_index_at_restart; ++fi.time_index) {
        ++fi.num_steps;
    }
}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::prepare_netcdf_file(
    ::kynema_sgf_tests::CountingTurbine& /*fi*/)
{}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::write_velocity_data(
    const ::kynema_sgf_tests::CountingTurbine& fi)
{
    fi.velocity_records.push_back(fi.time_index / fi.num_substeps);
}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::read_velocity_data(
    ::kynema_sgf_tests::CountingTurbine& /*fi*/,
    const ncutils::NCFile& /*ncf*/,
    const size_t /*tid*/)
{}

template class ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>;

} // namespace ext_turb

namespace kynema_sgf_tests {
namespace {

using CountingIface =
    ::ext_turb::ExtTurbIface<CountingTurbine, CountingSolverData>;
using CountingViewData =
    ExtTurbineViewData<CountingTurbine, CountingSolverData>;

class ExtTurbIfaceTest : public MeshTest
{
protected:
    //! Turbine with the hub as its only point, two solver steps per CFD step
    void setup_turbine(
        CountingViewData& data,
        CountingIface& iface,
        const ::ext_turb::SimMode sim_mode,
        const int cfd_step_at_restart)
    {
        init_hub_first_points(data, 0, 2, 0);
        data.info().is_root_proc = true;

        auto& fi = data.meta().ext_data;
        fi.tlabel = "T0";
        fi.tid_global = 0;
        fi.sim_mode = sim_mode;
        fi.dt_cfd = 0.5_rt;
        fi.dt_ext = 0.25_rt;
        fi.num_substeps = m_num_substeps;
        fi.chkpt_interval = 0;
        fi.time_index_at_restart = cfd_step_at_restart * m_num_substeps;

        data.meta().ext_ptr = &iface;
        iface.register_turbine(fi);
        iface.init_turbine(fi.tid_local);
    }

    //! What the actuator does once the flow is initialised: the loads of the
    //! turbine at the initial time (Actuator::post_init_actions)
    static void post_init(CountingViewData& data)
    {
        kynema_sgf::actuator::external::ext_loads(data);
    }

    //! What the actuator does before a CFD time step: the loads, then the
    //! advance of the turbine (Actuator::pre_advance_work)
    static void pre_advance(CountingViewData& data)
    {
        kynema_sgf::actuator::external::ext_loads(data);
        kynema_sgf::actuator::external::ext_step(data);
    }

    const int m_num_substeps{2};
};

} // namespace

TEST_F(ExtTurbIfaceTest, new_turbine_starts_with_initial_solution)
{
    initialize_mesh();
    CountingIface iface(sim());
    CountingViewData data;
    setup_turbine(data, iface, ::ext_turb::SimMode::init, 0);
    const auto& fi = data.meta().ext_data;

    post_init(data);
    EXPECT_EQ(fi.num_initial_solutions, 1);
    EXPECT_EQ(fi.num_steps, 0);
    EXPECT_EQ(fi.time_index, 0);
    EXPECT_TRUE(fi.velocity_records.empty());

    // The first two CFD time steps
    pre_advance(data);
    pre_advance(data);
    EXPECT_EQ(fi.num_initial_solutions, 1);
    EXPECT_EQ(fi.num_steps, 2 * m_num_substeps);
    EXPECT_EQ(fi.time_index, 2 * m_num_substeps);
    ASSERT_EQ(static_cast<int>(fi.velocity_records.size()), 2);
    EXPECT_EQ(fi.velocity_records[0], 0);
    EXPECT_EQ(fi.velocity_records[1], 1);
}

TEST_F(ExtTurbIfaceTest, restarted_turbine_has_no_second_initial_solution)
{
    const int cfd_step = 80;
    initialize_mesh();
    CountingIface iface(sim());
    CountingViewData data;
    setup_turbine(data, iface, ::ext_turb::SimMode::restart, cfd_step);
    const auto& fi = data.meta().ext_data;
    ASSERT_EQ(fi.time_index, cfd_step * m_num_substeps);

    // The checkpoint holds the turbine at the time of the flow: no initial
    // solution, which would start the extrapolation history of the solver
    // over
    post_init(data);
    EXPECT_EQ(fi.num_initial_solutions, 0);
    EXPECT_EQ(fi.num_steps, 0);
    EXPECT_EQ(fi.time_index, cfd_step * m_num_substeps);
    EXPECT_TRUE(fi.velocity_records.empty());

    // First CFD time step after the restart
    pre_advance(data);
    EXPECT_EQ(fi.num_initial_solutions, 0);
    EXPECT_EQ(fi.num_steps, m_num_substeps);
    EXPECT_EQ(fi.time_index, (cfd_step + 1) * m_num_substeps);
    ASSERT_EQ(static_cast<int>(fi.velocity_records.size()), 1);
    EXPECT_EQ(fi.velocity_records[0], cfd_step);
}

TEST_F(ExtTurbIfaceTest, replayed_turbine_stays_at_restart_time)
{
    const int cfd_step = 80;
    initialize_mesh();
    CountingIface iface(sim());
    CountingViewData data;
    setup_turbine(data, iface, ::ext_turb::SimMode::replay, cfd_step);
    const auto& fi = data.meta().ext_data;
    ASSERT_EQ(fi.time_index, cfd_step * m_num_substeps);
    ASSERT_EQ(fi.num_initial_solutions, 1);
    ASSERT_EQ(fi.num_steps, cfd_step * m_num_substeps);

    // The replay has brought the turbine to the time of the flow. One more
    // step here would leave it one CFD time step ahead for the whole run.
    post_init(data);
    EXPECT_EQ(fi.num_initial_solutions, 1);
    EXPECT_EQ(fi.num_steps, cfd_step * m_num_substeps);
    EXPECT_EQ(fi.time_index, cfd_step * m_num_substeps);
    EXPECT_TRUE(fi.velocity_records.empty());

    // First CFD time step after the restart
    pre_advance(data);
    EXPECT_EQ(fi.num_steps, (cfd_step + 1) * m_num_substeps);
    EXPECT_EQ(fi.time_index, (cfd_step + 1) * m_num_substeps);
    ASSERT_EQ(static_cast<int>(fi.velocity_records.size()), 1);
    EXPECT_EQ(fi.velocity_records[0], cfd_step);
}

} // namespace kynema_sgf_tests
