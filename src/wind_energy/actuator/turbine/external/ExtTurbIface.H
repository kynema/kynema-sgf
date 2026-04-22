#ifndef EXTTURBIFACE_H
#define EXTTURBIFACE_H

#include <map>
#include <vector>
#include "amr-wind/core/ExtSolver.H"
#include "amr-wind/wind_energy/actuator/turbine/external/external_base_types.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace ncutils {
class NCFile;
}

namespace amr_wind {
class CFDSim;
}

namespace ext_turb {

template <typename SolverTurbine, typename SolverData>
class ExtTurbIface
    : public ::amr_wind::ExtSolver::Register<
          ExtTurbIface<SolverTurbine, SolverData>>
{
public:
    static std::string identifier() { return ext_id<SolverData>(); }

    explicit ExtTurbIface(const ::amr_wind::CFDSim& sim);

    /* -- General for all external turbines -- */
    int register_turbine(SolverTurbine& data);

    void init_turbine(const int local_id);

    void advance_turbine(const int local_id);

    int num_local_turbines() const
    {
        return static_cast<int>(m_turbine_data.size());
    }

    /* -- Specific to each type of external turbine -- */

    ~ExtTurbIface() override;

    void parse_inputs(
        const amr_wind::CFDSim& /*sim*/, const std::string& /*inp_name*/);

    void init_solution(const int local_id);

    void get_hub_stats(const int local_id);

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
        const size_t tid);

    //! Global to local index lookup map
    std::map<int, int> m_turbine_map;

    std::vector<SolverTurbine*> m_turbine_data;

    amrex::Real m_dt_cfd{0.0_rt};
    amrex::Real m_start_time{0.0_rt};
    amrex::Real m_stop_time{-1.0_rt};
    ::ext_turb::SimMode m_sim_mode{::ext_turb::SimMode::init};

    SolverData m_solver_data;

    bool m_is_initialized{false};
};

// General implementations of some functions
template <typename SolverTurbine, typename SolverData>
ExtTurbIface<SolverTurbine, SolverData>::ExtTurbIface(
    const amr_wind::CFDSim& /*unused*/)
{}

template <typename SolverTurbine, typename SolverData>
int ExtTurbIface<SolverTurbine, SolverData>::register_turbine(
    SolverTurbine& data)
{
    BL_PROFILE("amr-wind::ExtTurbIface::register_turbine");
    AMREX_ALWAYS_ASSERT(!m_is_initialized);
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
    AMREX_ALWAYS_ASSERT(local_id < static_cast<int>(m_turbine_data.size()));
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
        break;
    }

    case ::ext_turb::SimMode::restart: {
        ext_restart_turbine(fi);
        prepare_netcdf_file(fi);
        break;
    }
    }
}

template <typename SolverTurbine, typename SolverData>
void ExtTurbIface<SolverTurbine, SolverData>::advance_turbine(
    const int local_id)
{
    BL_PROFILE("amr-wind::ExtTurbIface::advance_turbine");
    AMREX_ASSERT(local_id < static_cast<int>(m_turbine_data.size()));

    auto& fi = *m_turbine_data[local_id];
    AMREX_ASSERT(!fi.is_solution0);
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
                << std::endl << std::endl;
            // clang-format on
        }
    }

    write_velocity_data(fi);
    for (int i = 0; i < fi.num_substeps; ++i, ++fi.time_index) {
        do_turbine_step(fi);
    }

    if (fi.chkpt_interval > 0 &&
        (fi.time_index / fi.num_substeps) % fi.chkpt_interval == 0) {
        write_turbine_checkpoint(fi.tid_local);
    }
}

} // namespace ext_turb

#endif /* EXTTURBIFACE_H */
