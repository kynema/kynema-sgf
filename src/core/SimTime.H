#ifndef SIMTIME_H
#define SIMTIME_H

#include "AMReX_REAL.H"
#include "AMReX_GpuQualifiers.H"
#include "AMReX_Extension.H"
#include "AMReX_Vector.H"
#include <algorithm>
#include <cmath>

using namespace amrex::literals;

namespace {
AMREX_FORCE_INLINE amrex::Real get_enforced_dt_for_output(
    const amrex::Real dt,
    const amrex::Real cur_time,
    const amrex::Real interval,
    const amrex::Real tol)
{
    // Tolerance for being slightly under is relative to requested interval
    return amrex::min<amrex::Real>(
        dt, ((std::floor((cur_time / interval) + tol) + 1.0_rt) * interval) -
                cur_time);
}
} // namespace

namespace kynema_sgf {

/** Time manager for simulations
 *  \ingroup core
 */
class SimTime
{
public:
    static constexpr int max_time_states = 3;

    SimTime() = default;

    ~SimTime() = default;

    /** Advance to a new time state and update internal time states
     *
     *  Return a boolean flag indicating if the simulation should continue.
     */
    bool new_timestep();

    /** Return true if simulation should continue
     */
    [[nodiscard]] bool continue_simulation() const;

    /** Return true if mesh refinement should be performed at this timestep
     */
    [[nodiscard]] bool do_regrid() const;

    /** Return true if plot files should be written at this timestep
     */
    [[nodiscard]] bool write_plot_file() const;

    /** Return true if checkpoint files should be written at this timestep
     */
    [[nodiscard]] bool write_checkpoint() const;

    /** Return true if plot files should be written at this timestep
     */
    [[nodiscard]] bool write_last_plot_file() const;

    /** Return true if checkpoint files should be written at this timestep
     */
    [[nodiscard]] bool write_last_checkpoint() const;

    /** Set current CFL and update timestep based on CFL components
     *
     */
    void set_current_cfl(
        amrex::Real conv_cfl, amrex::Real diff_cfl, amrex::Real src_cfl);

    /** Use results of CFL and timestepping parameters to advance time to new,
     * output time information
     */
    void advance_time();

    /** Set start time and index based on checkpoint/restart file
     *
     *  \param tidx Time index from checkpoint file
     *  \param time Time read from checkpoint file
     */
    void set_restart_time(int tidx, amrex::Real time);

    AMREX_FORCE_INLINE
    amrex::Real delta_t() const { return m_dt[0]; }

    AMREX_FORCE_INLINE
    amrex::Real& delta_t() { return m_dt[0]; }

    AMREX_FORCE_INLINE
    amrex::Real delta_t_nm1() const { return m_dt[1]; }

    AMREX_FORCE_INLINE
    amrex::Real delta_t_nm2() const { return m_dt[2]; }

    AMREX_FORCE_INLINE
    amrex::Real& delta_t_nm1() { return m_dt[1]; }

    AMREX_FORCE_INLINE
    amrex::Real& delta_t_nm2() { return m_dt[2]; }

    AMREX_FORCE_INLINE
    amrex::Real current_time() const { return m_cur_time; }

    AMREX_FORCE_INLINE
    amrex::Real& current_time() { return m_cur_time; }

    AMREX_FORCE_INLINE
    amrex::Real new_time() const { return m_new_time; }

    AMREX_FORCE_INLINE
    amrex::Real max_cfl() const { return m_max_cfl; }

    AMREX_FORCE_INLINE
    amrex::Real get_minimum_enforce_dt_abs_tol() const
    {
        return m_force_dt_abs_tol;
    }

    AMREX_FORCE_INLINE
    int time_index() const { return m_time_index; }

    AMREX_FORCE_INLINE
    int& time_index() { return m_time_index; }

    AMREX_FORCE_INLINE
    bool adaptive_timestep() const { return m_adaptive; }

    AMREX_FORCE_INLINE
    bool use_force_cfl() const { return m_use_force_cfl; }

    AMREX_FORCE_INLINE
    int regrid_interval() const { return m_regrid_interval; }

    AMREX_FORCE_INLINE
    amrex::Real start_time() const { return m_start_time; }

    AMREX_FORCE_INLINE
    amrex::Real stop_time() const { return m_stop_time; }

    AMREX_FORCE_INLINE
    int start_time_index() const { return m_start_time_index; }

    AMREX_FORCE_INLINE
    int stop_time_index() const { return m_stop_time_index; }

    AMREX_FORCE_INLINE
    int chkpt_interval() const { return m_chkpt_interval; }

    //! Read user defined options from input file
    void parse_parameters();

    //! Get timestep considerations from postprocessing instances
    void add_postproc_dt_parameters(
        const bool enforce_dt_flag,
        const amrex::Real enforce_dt_tol,
        const amrex::Real time_interval,
        const amrex::Real time_delay)
    {
        m_postprocess_enforce_dt.emplace_back(enforce_dt_flag);
        m_postprocess_enforce_dt_tol.emplace_back(enforce_dt_tol);
        m_postprocess_time_interval.emplace_back(time_interval);
        m_postprocess_time_delay.emplace_back(time_delay);
    }

    //! Calculate minimum enforce dt tolerance based on all considerations
    void calculate_minimum_enforce_dt_abs_tol();

    //! Deactivate ending parameters (overset should be externally controlled)
    void override_simulation_end_parameters()
    {
        m_stop_time_index = -1;
        m_stop_time = -1.0_rt;
    }

    /** Return true if TinyProfiler output should be sent to the log file at
     * this timestep
     */
    [[nodiscard]] bool output_profiling_info() const;

private:
    //! Timestep sizes
    amrex::Vector<amrex::Real> m_dt =
        amrex::Vector<amrex::Real>(max_time_states, 0.0_rt);

    //! Current simulation time
    amrex::Real m_cur_time{0.0_rt};

    //! Time where a new timestep will advance to
    amrex::Real m_new_time{0.0_rt};

    //! Starting simulation time (to track restart)
    amrex::Real m_start_time{0.0_rt};

    //! Max time for simulation
    amrex::Real m_stop_time{-1.0_rt};

    //! Delayed adaptive stepping
    amrex::Real m_delay_time{-1.0_rt};

    //! Maximum CFL constraint
    amrex::Real m_max_cfl{1.0_rt};

    //! Fixed delta_t indicated by user
    amrex::Real m_fixed_dt{-1.0_rt};

    //! Initial delta_t indicated by user
    amrex::Real m_initial_dt{-1.0_rt};

    //! Maximum delta_t indicated by user
    amrex::Real m_max_dt{-1.0_rt};

    //! Minimum delta_t indicated by user
    amrex::Real m_min_dt{-1.0_rt};

    //! Current CFL
    amrex::Real m_current_cfl{0.0_rt};

    //! Current CFL categories
    amrex::Real m_conv_cfl{0.0_rt};
    amrex::Real m_diff_cfl{0.0_rt};
    amrex::Real m_src_cfl{0.0_rt};

    //! delta_t calculated using CFL constraint
    amrex::Real m_dt_calc{0.0_rt};

    //! Initial reduction in timestep size for startup
    amrex::Real m_init_shrink{0.1_rt};

    //! Maximum growth of dt between timesteps
    amrex::Real m_dt_growth{0.1_rt};

    //! Time interval for writing checkpoint/restart files
    amrex::Real m_chkpt_t_interval{-1.0_rt};

    //! Time delay for checkpoint/restart output
    amrex::Real m_chkpt_t_delay{0.0_rt};

    //! Relative (to dt) tolerance for checkpoint time interval output
    amrex::Real m_chkpt_t_tol{1.0e-2_rt};

    //! Relative (to chkpt_t_interval) tolerance for enforcing dt
    amrex::Real m_force_chkpt_tol{1.0e-3_rt};

    //! Time interval for plot file output
    amrex::Real m_plt_t_interval{-1.0_rt};

    //! Time delay for plot file output
    amrex::Real m_plt_t_delay{0.0_rt};

    //! Relative (to dt) tolerance for plot time interval output
    amrex::Real m_plt_t_tol{1.0e-2_rt};

    //! Relative (to plt_t_interval) tolerance for enforcing dt
    amrex::Real m_force_plt_tol{1.0e-3_rt};

    //! Smallest absolute tolerance for enforcing dt, considering all
    //! tolerance criteria. This is simply a high initial value; the code
    //! calculates the actual tolerance with min() operations
    amrex::Real m_force_dt_abs_tol{1.0e8_rt};

    //! Initial starting time for generating checkpoint files (time-based)
    amrex::Real m_chkpt_start_time{0.0_rt};

    //! Initial starting time for generating output files (time-based)
    amrex::Real m_plt_start_time{0.0_rt};

    //! Counter for the number of timesteps since start of simulation
    int m_time_index{0};

    //! Initial starting time index
    int m_start_time_index{0};

    //! Initial starting time index for generating checkpoint files
    int m_chkpt_start_index{0};

    //! Initial starting time index for generating output files
    int m_plt_start_index{0};

    //! Initial starting time index for regriding
    int m_regrid_start_index{0};

    //! Maximum timesteps for simulation
    int m_stop_time_index{-1};

    //! Time step interval for writing checkpoint/restart files
    int m_chkpt_interval{-1};

    //! Time step delay for checkpoint/restart output
    //! (this means that output only happens after the delay)
    int m_chkpt_delay{0};

    //! Time step interval for plot file output
    int m_plt_interval{-1};

    //! Time step delay for plot file output
    int m_plt_delay{0};

    //! Time interval for regridding
    int m_regrid_interval{-1};

    //! Time step interval for profiling output
    int m_profiling_interval{-1};

    //! Verbosity
    int m_verbose{0};

    //! Flag indicating whether timestep is adaptive
    bool m_adaptive{false};

    //! Flag indicating if this is initialization
    bool m_is_init{true};

    //! Flag indicating if forcing should be included in CFL calculation
    bool m_use_force_cfl{true};

    //! Bool for if checkpoint time interval should be forced
    bool m_force_chkpt_dt{false};

    //! Bool for if plt time interval should be forced
    bool m_force_plt_dt{false};

    //! Timestep considerations from postprocessing instances
    std::vector<bool> m_postprocess_enforce_dt;
    amrex::Vector<amrex::Real> m_postprocess_enforce_dt_tol;
    amrex::Vector<amrex::Real> m_postprocess_time_interval;
    amrex::Vector<amrex::Real> m_postprocess_time_delay;
};

} // namespace kynema_sgf

#endif /* SIMTIME_H */
