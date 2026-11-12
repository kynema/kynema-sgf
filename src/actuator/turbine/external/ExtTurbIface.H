#ifndef EXTTURBIFACE_H
#define EXTTURBIFACE_H

#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "src/CFDSim.H"
#include "src/core/ExtSolver.H"
#include "src/actuator/turbine/external/external_base_types.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace ncutils {
class NCFile;
}

namespace ext_turb {

/** Number of time steps that make up a time interval, which has to be whole
 *
 *  Gives the steps of the external solver in one CFD time step and the CFD
 *  time steps up to the start time of a replay. The ratio is computed in
 *  floating point, where 0.3 / 0.1 is 2.9999999999999996: the nearest integer
 *  is taken, and the run stops when the ratio is further from it than
 *  rounding explains.
 *
 *  \param interval Length of the interval
 *  \param dt Time step
 *  \param interval_name What the interval is, the start of the abort message
 *  \param dt_name What the time steps are, in the plural
 */
inline int num_whole_steps(
    const amrex::Real interval,
    const amrex::Real dt,
    const std::string& interval_name,
    const std::string& dt_name)
{
    // The rounding error of the ratio is a few units of machine precision
    constexpr amrex::Real tol =
        1.0e4_rt * std::numeric_limits<amrex::Real>::epsilon();
    constexpr auto max_steps =
        static_cast<amrex::Real>(std::numeric_limits<int>::max());

    const amrex::Real ratio = (dt > 0.0_rt) ? interval / dt : -1.0_rt;
    const amrex::Real nearest = std::round(ratio);
    if (!(ratio >= 0.0_rt) || !(nearest < max_steps) ||
        (std::abs(ratio - nearest) > tol * ratio)) {
        std::ostringstream msg;
        msg << std::setprecision(12) << interval_name << " (" << interval
            << ") is not a whole number of " << dt_name << " (" << dt << ")";
        amrex::Abort(msg.str());
    }
    return static_cast<int>(nearest);
}

template <typename SolverTurbine, typename SolverData>
class ExtTurbIface
    : public ::kynema_sgf::ExtSolver::Register<
          ExtTurbIface<SolverTurbine, SolverData>>
{
public:
    static std::string identifier() { return ext_id<SolverData>(); }

    explicit ExtTurbIface(const ::kynema_sgf::CFDSim& sim);

    /* -- General for all external turbines -- */
    int register_turbine(SolverTurbine& data);

    void init_turbine(int local_id);

    void advance_turbine(int local_id);

    void advance_aero_loads(int local_id);

    [[nodiscard]] int num_local_turbines() const
    {
        return static_cast<int>(m_turbine_data.size());
    }

    /* -- Specific to each type of external turbine -- */

    ~ExtTurbIface() override;

    void parse_inputs(
        const kynema_sgf::CFDSim& /*sim*/, const std::string& /*inp_name*/);

    void init_solution(int local_id);

    void get_hub_stats(int local_id);

    void calc_aero_loads(SolverTurbine& fi);

    void do_turbine_step(SolverTurbine& fi);

    void write_turbine_checkpoint(int& tid);

protected:
    void allocate_ext_turbines();

    void ext_init_turbine(SolverTurbine& /*fi*/);

    void ext_restart_turbine(SolverTurbine& /*fi*/);

    void ext_replay_turbine(SolverTurbine& /*fi*/);

    void prepare_netcdf_file(SolverTurbine& /*unused*/);

    void write_velocity_data(const SolverTurbine& /*unused*/);

    static void read_velocity_data(
        SolverTurbine& /*unused*/,
        const ncutils::NCFile& /*unused*/,
        size_t tid);

    //! Global to local index lookup map
    std::map<int, int> m_turbine_map;

    std::vector<SolverTurbine*> m_turbine_data;

    amrex::Real m_dt_cfd{0.0_rt};
    amrex::Real m_start_time{0.0_rt};
    amrex::Real m_stop_time{-1.0_rt};
    ::ext_turb::SimMode m_sim_mode{::ext_turb::SimMode::init};

    SolverData m_solver_data;

    bool m_is_initialized{false};

    //! Time of the CFD solver, which says when a checkpoint is written
    const ::kynema_sgf::SimTime& m_time;
};

// General implementations of some functions
template <typename SolverTurbine, typename SolverData>
ExtTurbIface<SolverTurbine, SolverData>::ExtTurbIface(
    const kynema_sgf::CFDSim& sim)
    : m_time(sim.time())
{}

template <typename SolverTurbine, typename SolverData>
int ExtTurbIface<SolverTurbine, SolverData>::register_turbine(
    SolverTurbine& data)
{
    BL_PROFILE("kynema-sgf::ExtTurbIface::register_turbine");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !m_is_initialized,
        "m_is_initialized = " + std::to_string(m_is_initialized));
    const int local_id = static_cast<int>(m_turbine_data.size());
    const int gid = data.tid_global;
    m_turbine_map[gid] = local_id;
    data.tid_local = local_id;
    m_turbine_data.emplace_back(&data);

    return local_id;
}

template <typename SolverTurbine, typename SolverData>
void ExtTurbIface<SolverTurbine, SolverData>::init_turbine(const int local_id)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        local_id < static_cast<int>(m_turbine_data.size()),
        "local_id = " + std::to_string(local_id) +
            ", m_turbine_data.size() = " +
            std::to_string(m_turbine_data.size()));
    if (!m_is_initialized) {
        allocate_ext_turbines();
    }
    auto& fi = *m_turbine_data[local_id];

    switch (fi.sim_mode) {
    case ::ext_turb::SimMode::init: {
        ext_init_turbine(fi);
        prepare_netcdf_file(fi);
        break;
    }

    case ::ext_turb::SimMode::replay: {
        ext_init_turbine(fi);
        ext_replay_turbine(fi);
        fi.is_solution0 = false;
        break;
    }

    case ::ext_turb::SimMode::restart: {
        ext_restart_turbine(fi);
        prepare_netcdf_file(fi);
        fi.is_solution0 = false;
        break;
    }
    }
}

template <typename SolverTurbine, typename SolverData>
void ExtTurbIface<SolverTurbine, SolverData>::advance_aero_loads(
    const int local_id)
{
    BL_PROFILE("kynema-sgf::ExtTurbIface::advance_aero_loads");
    AMREX_ASSERT(local_id < static_cast<int>(m_turbine_data.size()));

    auto& fi = *m_turbine_data[local_id];
    calc_aero_loads(fi);
}

template <typename SolverTurbine, typename SolverData>
void ExtTurbIface<SolverTurbine, SolverData>::advance_turbine(
    const int local_id)
{
    BL_PROFILE("kynema-sgf::ExtTurbIface::advance_turbine");
    AMREX_ASSERT(local_id < static_cast<int>(m_turbine_data.size()));

    auto& fi = *m_turbine_data[local_id];
    AMREX_ASSERT(!fi.is_solution0);
    write_velocity_data(fi);
    {
        // Default is off, unless a turbine model reads, populates stop_time
        const auto& tmax = fi.stop_time;
        const auto& telapsed = (fi.time_index + fi.num_substeps) * fi.dt_ext;
        if (telapsed > (tmax + std::numeric_limits<float>::epsilon()) &&
            fi.stop_time > 0.0_rt) {
            // clang-format off
            amrex::OutStream()
                << "\nWARNING: ExtTurbIface:\n"
                << "  Elapsed simulation time will exceed max "
                << "time set for External Turbine Solver"
                << "\n\n";
            // clang-format on
        }
    }

    for (int i = 0; i < fi.num_substeps; ++i, ++fi.time_index) {
        do_turbine_step(fi);
    }

    // The checkpoint of the external solver goes with the one the CFD solver
    // writes at the end of this time step: on its step or time interval,
    // after its delay, counted from the restart, and at the last step of the
    // run when that one is off the interval
    if (m_time.write_checkpoint() ||
        (!m_time.continue_simulation() && m_time.write_last_checkpoint())) {
        write_turbine_checkpoint(fi.tid_local);
    }
}

} // namespace ext_turb

#endif /* EXTTURBIFACE_H */
