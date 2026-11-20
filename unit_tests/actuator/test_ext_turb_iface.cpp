#include <array>
#include <memory>
#include <string>
#include <vector>
#include "ks_test_utils/MeshTest.H"
#include "test_ext_turbine_utils.H"
#include "src/actuator/turbine/external/ExtTurbIface.H"
#include "src/actuator/turbine/external/turbine_external_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

/** Turbine of a stand-in external solver without physics
 *
 *  Initial solutions, steps, checkpoints and velocity records are counted.
 */
struct CountingTurbine : public ::ext_turb::ExternalTurbine
{
    float* position_at_vel(int /*dir*/) const override { return nullptr; }
    float* solid_velocity(int /*dir*/) const override { return nullptr; }
    float* fluid_velocity(int dir) const override { return &fluid_vel[dir]; }
    float* force(int dir) const override { return &hub_force[dir]; }
    float* position_at_force(int /*dir*/) const override { return nullptr; }
    float* chord_at_force() const override { return chord.data(); }
    float* orientation() const override { return nullptr; }

    int length_fluid_velocity(int /*dir*/) const override { return num_pts(); }
    int length_force(int /*dir*/) const override { return num_pts(); }
    int length_position_at_force(int /*dir*/) const override
    {
        return num_pts();
    }
    int length_orientation() const override { return 0; }
    int num_vel_pts_blade() const override { return model_num_pts_blade; }

    //! Hub and blade points of the solver model
    int num_pts() const { return 1 + (model_num_blades * model_num_pts_blade); }

    //! Blades of the solver model and points per blade; default: hub only
    int model_num_blades{0};
    int model_num_pts_blade{0};

    //! Chord at the points of the model, for the tests that set blades
    mutable std::vector<float> chord;

    //! Solver step count at the CFD restart
    int time_index_at_restart{0};

    int num_initial_solutions{0};
    int num_steps{0};
    int num_checkpoints{0};

    //! CFD time step of each velocity record written
    mutable std::vector<int> velocity_records;

    //! CFD time step of each checkpoint written
    std::vector<int> checkpoint_steps;

    //! Hub velocity and force
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
    ::kynema_sgf_tests::CountingSolverData>::
    parse_inputs(
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
    ::kynema_sgf_tests::CountingSolverData>::
    get_hub_stats(const int /*local_id*/)
{}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::
    calc_aero_loads(::kynema_sgf_tests::CountingTurbine& /*fi*/)
{}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::
    do_turbine_step(::kynema_sgf_tests::CountingTurbine& fi)
{
    ++fi.num_steps;
}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::write_turbine_checkpoint(int& tid)
{
    auto& fi = *m_turbine_data[tid];
    ++fi.num_checkpoints;
    fi.checkpoint_steps.push_back(fi.time_index / fi.num_substeps);
}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::
    ext_init_turbine(::kynema_sgf_tests::CountingTurbine& fi)
{
    fi.time_index = 0;
    // Number of blades from the solver model, as in OpenFAST
    fi.num_blades = fi.model_num_blades;
}

// Checkpoint read: solver at the step of the checkpoint
template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::
    ext_restart_turbine(::kynema_sgf_tests::CountingTurbine& fi)
{
    fi.time_index = fi.time_index_at_restart;
}

// Replay: initial solution, then the steps up to the restart time
template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::
    ext_replay_turbine(::kynema_sgf_tests::CountingTurbine& fi)
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
    ::kynema_sgf_tests::CountingSolverData>::
    prepare_netcdf_file(::kynema_sgf_tests::CountingTurbine& /*fi*/)
{}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::
    write_velocity_data(const ::kynema_sgf_tests::CountingTurbine& fi)
{
    fi.velocity_records.push_back(fi.time_index / fi.num_substeps);
}

template <>
void ExtTurbIface<
    ::kynema_sgf_tests::CountingTurbine,
    ::kynema_sgf_tests::CountingSolverData>::
    read_velocity_data(
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
    //! Interface of a test, owned by the fixture: the view data point to it
    CountingIface& make_iface()
    {
        m_iface = std::make_unique<CountingIface>(sim());
        return *m_iface;
    }

    //! View data of a test, owned by the fixture: the interface points to it
    CountingViewData& make_data()
    {
        m_data = std::make_unique<CountingViewData>();
        return *m_data;
    }

    //! Hub-only turbine, two solver steps per CFD step
    void setup_turbine(
        CountingViewData& data,
        CountingIface& iface,
        const ::ext_turb::SimMode sim_mode,
        const int cfd_step_at_restart) const
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
        fi.time_index_at_restart = cfd_step_at_restart * m_num_substeps;

        data.meta().ext_ptr = &iface;
        iface.register_turbine(fi);
        iface.init_turbine(fi.tid_local);
    }

    //! Actuator::post_init_actions: loads at the initial time
    static void post_init(CountingViewData& data)
    {
        kynema_sgf::actuator::external::ext_loads(data);
    }

    //! Actuator::pre_advance_work: loads, then advance
    static void pre_advance(CountingViewData& data)
    {
        kynema_sgf::actuator::external::ext_loads(data);
        kynema_sgf::actuator::external::ext_step(data);
    }

    //! Run inputs: max_step, fixed dt of the turbine. Tests add checkpoint
    //! inputs.
    void time_inputs(const int max_step)
    {
        populate_parameters();
        amrex::ParmParse pp("time");
        pp.add("stop_time", 100.0_rt);
        pp.add("fixed_dt", 0.5_rt);
        pp.add("max_step", max_step);
    }

    //! IOManager::write_checkpoint_file, as far as the turbine goes
    void write_cfd_checkpoint(
        CountingViewData& data, std::vector<int>& cfd_checkpoint_steps)
    {
        cfd_checkpoint_steps.push_back(sim().time().time_index());
        kynema_sgf::actuator::ops::external_write_solver_checkpoint(data);
    }

    /** incflo::Evolve, as far as checkpoints go
     *
     *  Advance after the new time is set; checkpoint at the end of a step and
     *  after the loop (write_final_checkpoint). stop_step > 0: stop requested
     *  at the end of that step (as RANSConvergence). Returns the CFD checkpoint
     *  steps.
     */
    std::vector<int>
    run_time_loop(CountingViewData& data, const int stop_step = 0)
    {
        // Loads at the initial time, or nothing after a restart
        post_init(data);

        std::vector<int> cfd_checkpoint_steps;
        auto& time = sim().time();
        while (time.new_timestep()) {
            time.set_current_cfl(0.25_rt, 0.0_rt, 0.0_rt);
            time.advance_time();
            pre_advance(data);
            if (time.time_index() == stop_step) {
                time.request_stop("test");
            }
            if (time.write_checkpoint()) {
                write_cfd_checkpoint(data, cfd_checkpoint_steps);
            }
        }
        if (time.write_final_checkpoint()) {
            write_cfd_checkpoint(data, cfd_checkpoint_steps);
        }
        return cfd_checkpoint_steps;
    }

    const int m_num_substeps{2};

    std::unique_ptr<CountingIface> m_iface;
    std::unique_ptr<CountingViewData> m_data;
};

} // namespace

TEST_F(ExtTurbIfaceTest, new_turbine_starts_with_initial_solution)
{
    initialize_mesh();
    auto& iface = make_iface();
    auto& data = make_data();
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

TEST_F(ExtTurbIfaceTest, fllc_data_follow_blade_count_of_solver)
{
    const int num_blades_input = 2;
    const int num_blades_model = 3;
    const int num_pts_blade = 4;
    initialize_mesh();
    auto& iface = make_iface();
    auto& data = make_data();
    data.info().is_root_proc = true;
    data.info().root_proc = amrex::ParallelDescriptor::MyProc();

    // From the inputs: num_blades = 2, FLLC data for two blades
    auto& meta = data.meta();
    meta.num_blades = num_blades_input;
    meta.num_pts_blade = num_pts_blade;
    meta.num_pts_tower = 0;
    meta.eps_inp = kynema_sgf::vs::Vector::one();
    meta.fllc.resize(num_blades_input);

    // The model of the solver has three blades
    auto& fi = meta.ext_data;
    fi.tlabel = "T0";
    fi.tid_global = 0;
    fi.sim_mode = ::ext_turb::SimMode::init;
    fi.num_pts_tower = 0;
    fi.model_num_blades = num_blades_model;
    fi.model_num_pts_blade = num_pts_blade;
    fi.chord.assign(fi.num_pts(), 1.0F);
    meta.ext_ptr = &iface;
    iface.register_turbine(fi);

    kynema_sgf::actuator::ops::init_data_op(data);

    // Blades of the model, each with its FLLC data
    EXPECT_EQ(meta.num_blades, num_blades_model);
    EXPECT_EQ(static_cast<int>(meta.blades.size()), num_blades_model);
    EXPECT_EQ(static_cast<int>(meta.fllc.size()), num_blades_model);
}

TEST_F(ExtTurbIfaceTest, restarted_turbine_has_no_second_initial_solution)
{
    const int cfd_step = 80;
    initialize_mesh();
    auto& iface = make_iface();
    auto& data = make_data();
    setup_turbine(data, iface, ::ext_turb::SimMode::restart, cfd_step);
    const auto& fi = data.meta().ext_data;
    ASSERT_EQ(fi.time_index, cfd_step * m_num_substeps);

    // Checkpoint at the flow time: no initial solution (it would reset the
    // solver extrapolation history)
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
    auto& iface = make_iface();
    auto& data = make_data();
    setup_turbine(data, iface, ::ext_turb::SimMode::replay, cfd_step);
    const auto& fi = data.meta().ext_data;
    ASSERT_EQ(fi.time_index, cfd_step * m_num_substeps);
    ASSERT_EQ(fi.num_initial_solutions, 1);
    ASSERT_EQ(fi.num_steps, cfd_step * m_num_substeps);

    // Replay ends at the flow time: no step here
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

// A CFD checkpoint needs a solver checkpoint of the same step: the tests
// below compare the two

TEST_F(ExtTurbIfaceTest, checkpoint_at_last_step_off_the_interval)
{
    time_inputs(7);
    {
        amrex::ParmParse pp("time");
        pp.add("checkpoint_interval", 5);
    }
    initialize_mesh();
    auto& iface = make_iface();
    auto& data = make_data();
    setup_turbine(data, iface, ::ext_turb::SimMode::init, 0);

    // On the interval at step 5, and at the end of the run
    const std::vector<int> expected{5, 7};
    EXPECT_EQ(run_time_loop(data), expected);
    EXPECT_EQ(data.meta().ext_data.checkpoint_steps, expected);
}

TEST_F(ExtTurbIfaceTest, checkpoint_on_time_interval)
{
    time_inputs(7);
    {
        amrex::ParmParse pp("time");
        pp.add("checkpoint_interval", -1);
        pp.add("checkpoint_time_interval", 1.0_rt);
    }
    initialize_mesh();
    auto& iface = make_iface();
    auto& data = make_data();
    setup_turbine(data, iface, ::ext_turb::SimMode::init, 0);

    // Every 1.0 of time with a time step of 0.5, and at the end of the run
    const std::vector<int> expected{2, 4, 6, 7};
    EXPECT_EQ(run_time_loop(data), expected);
    EXPECT_EQ(data.meta().ext_data.checkpoint_steps, expected);
}

TEST_F(ExtTurbIfaceTest, checkpoint_interval_counted_from_restart)
{
    const int cfd_step = 3;
    time_inputs(10);
    {
        amrex::ParmParse pp("time");
        pp.add("checkpoint_interval", 5);
    }
    initialize_mesh();
    // CFD time after reading the checkpoint of step 3
    sim().time().set_restart_time(cfd_step, 1.5_rt);
    auto& iface = make_iface();
    auto& data = make_data();
    setup_turbine(data, iface, ::ext_turb::SimMode::restart, cfd_step);

    // Five steps after the restart, and at the end of the run
    const std::vector<int> expected{8, 10};
    EXPECT_EQ(run_time_loop(data), expected);
    EXPECT_EQ(data.meta().ext_data.checkpoint_steps, expected);
}

TEST_F(ExtTurbIfaceTest, checkpoint_after_delay)
{
    time_inputs(12);
    {
        amrex::ParmParse pp("time");
        pp.add("checkpoint_interval", 5);
        pp.add("checkpoint_delay", 6);
    }
    initialize_mesh();
    auto& iface = make_iface();
    auto& data = make_data();
    setup_turbine(data, iface, ::ext_turb::SimMode::init, 0);

    // Step 5 is before the delay
    const std::vector<int> expected{10, 12};
    EXPECT_EQ(run_time_loop(data), expected);
    EXPECT_EQ(data.meta().ext_data.checkpoint_steps, expected);
}

TEST_F(ExtTurbIfaceTest, checkpoint_with_final_checkpoint_after_stop_request)
{
    time_inputs(10);
    {
        amrex::ParmParse pp("time");
        pp.add("checkpoint_interval", 5);
    }
    initialize_mesh();
    auto& iface = make_iface();
    auto& data = make_data();
    setup_turbine(data, iface, ::ext_turb::SimMode::init, 0);

    // Stop requested after the advance of step 3: final checkpoint, off the
    // interval, with the turbine
    const std::vector<int> expected{3};
    EXPECT_EQ(run_time_loop(data, 3), expected);
    EXPECT_EQ(data.meta().ext_data.checkpoint_steps, expected);
}

TEST_F(ExtTurbIfaceTest, advance_writes_no_checkpoint)
{
    time_inputs(5);
    {
        amrex::ParmParse pp("time");
        pp.add("checkpoint_interval", 1);
    }
    initialize_mesh();
    auto& iface = make_iface();
    auto& data = make_data();
    setup_turbine(data, iface, ::ext_turb::SimMode::init, 0);

    // Every step is a checkpoint step, but the advance (turbine root only)
    // writes none: only a CFD checkpoint does
    post_init(data);
    auto& time = sim().time();
    while (time.new_timestep()) {
        time.set_current_cfl(0.25_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        pre_advance(data);
    }
    EXPECT_EQ(data.meta().ext_data.num_steps, 5 * m_num_substeps);
    EXPECT_TRUE(data.meta().ext_data.checkpoint_steps.empty());
}

TEST_F(ExtTurbIfaceTest, checkpoint_only_on_root_process)
{
    time_inputs(5);
    initialize_mesh();
    auto& iface = make_iface();
    auto& data = make_data();
    setup_turbine(data, iface, ::ext_turb::SimMode::init, 0);

    data.info().is_root_proc = false;
    kynema_sgf::actuator::ops::external_write_solver_checkpoint(data);
    EXPECT_TRUE(data.meta().ext_data.checkpoint_steps.empty());

    data.info().is_root_proc = true;
    kynema_sgf::actuator::ops::external_write_solver_checkpoint(data);
    EXPECT_EQ(data.meta().ext_data.checkpoint_steps.size(), 1U);
}

} // namespace kynema_sgf_tests
