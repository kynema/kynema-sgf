#include "src/core/SimTime.H"

#include "AMReX_ParmParse.H"
#include "AMReX_Print.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

void SimTime::parse_parameters()
{
    // Initialize delta_t to negative values
    for (amrex::Real& i : m_dt) {
        i = -1.0_rt;
    }

    // Parse options in "time" namespace
    amrex::ParmParse pp("time");
    pp.query("stop_time", m_stop_time);
    pp.query("max_step", m_stop_time_index);
    pp.query("delay_time", m_delay_time);
    pp.query("fixed_dt", m_fixed_dt);
    pp.query("initial_dt", m_initial_dt);
    pp.query("max_dt", m_max_dt);
    pp.query("min_dt", m_min_dt);
    pp.query("init_shrink", m_init_shrink);
    pp.query("max_dt_growth", m_dt_growth);
    pp.query("cfl", m_max_cfl);
    pp.query("verbose", m_verbose);
    pp.query("regrid_interval", m_regrid_interval);
    pp.query("plot_interval", m_plt_interval);
    pp.query("plot_time_interval", m_plt_t_interval);
    pp.query("plot_delay", m_plt_delay);
    pp.query("plot_time_delay", m_plt_t_delay);
    pp.query("enforce_plot_time_dt", m_force_plt_dt);
    pp.query("checkpoint_interval", m_chkpt_interval);
    pp.query("checkpoint_time_interval", m_chkpt_t_interval);
    pp.query("checkpoint_delay", m_chkpt_delay);
    pp.query("checkpoint_time_delay", m_chkpt_t_delay);
    pp.query("enforce_checkpoint_time_dt", m_force_chkpt_dt);
    pp.query("regrid_start", m_regrid_start_index);
    pp.query("plot_start", m_plt_start_index);
    pp.query("checkpoint_start", m_chkpt_start_index);
    pp.query("plot_start_time", m_plt_start_time);
    pp.query("checkpoint_start_time", m_chkpt_start_time);
    pp.query("use_force_cfl", m_use_force_cfl);
    pp.query("profiling_interval", m_profiling_interval);

    // Tolerances
    pp.query("plot_time_interval_reltol", m_plt_t_tol);
    pp.query("enforce_plot_dt_reltol", m_force_plt_tol);
    pp.query("checkpoint_time_interval_reltol", m_chkpt_t_tol);
    pp.query("enforce_checkpoint_dt_reltol", m_force_chkpt_tol);

    if (m_fixed_dt > 0.0_rt) {
        m_dt[0] = m_fixed_dt;
    } else {
        m_adaptive = true;
    }

    if (m_plt_interval > 0 && m_plt_t_interval > 0.0_rt) {
        amrex::Abort(
            "plot_interval and plot_time_interval are both specified. "
            "timestep- and time-based plotting should not be used together; "
            "please only specify one.");
    }

    if (m_chkpt_interval > 0 && m_chkpt_t_interval > 0.0_rt) {
        amrex::Abort(
            "checkpoint_interval and checkpoint_time_interval are both "
            "specified. timestep- and time-based checkpointing should not be "
            "used together; please only specify one.");
    }

    if (m_plt_t_interval <= 0.0_rt && m_force_plt_dt) {
        amrex::Abort(
            "enforce_plot_time_dt is true, but no plot time interval has been "
            "provided.");
    }

    if (m_chkpt_t_interval <= 0.0_rt && m_force_chkpt_dt) {
        amrex::Abort(
            "enforce_checkpoint_time_dt is true, but no checkpoint time "
            "interval has been provided.");
    }

    if (!m_adaptive && (m_force_plt_dt || m_force_chkpt_dt)) {
        amrex::Abort(
            "an output time interval has been specified to be enforced upon "
            "dt, but dt is not adaptive.");
    }

    if (m_force_plt_dt && m_force_chkpt_dt) {
        amrex::Print()
            << "WARNING: Time intervals will be enforced upon dt for both "
               "plotfiles and checkpoint files.";
        amrex::Print()
            << " -- If these time intervals are different and not factors of "
               "each other, inefficient behavior may result, depending on "
               "tolerances: dt may be shortened often and outputs may occur in "
               "consecutive timesteps.";
    }

    calculate_minimum_enforce_dt_abs_tol();
}

bool SimTime::new_timestep()
{
    bool continue_sim = continue_simulation();

    if (m_is_init && (m_verbose >= 0)) {
        amrex::Print() << "\nBegin simulation: " << '\n';
        if ((m_stop_time > 0) && (m_stop_time_index >= 0)) {
            amrex::Print() << "  Run until " << m_stop_time << " sec. or "
                           << m_stop_time_index << " timesteps" << '\n';
        } else if (m_stop_time > 0) {
            amrex::Print() << "  Run till " << m_stop_time << " seconds "
                           << '\n';
        } else if (m_stop_time_index >= 0) {
            amrex::Print() << "  Run for " << m_stop_time_index << " timesteps"
                           << '\n';
        }
        if (m_adaptive) {
            amrex::Print() << "  Adaptive timestepping with max. CFL = "
                           << m_max_cfl << '\n';
        } else {
            amrex::Print() << "  Fixed timestepping with dt = " << m_fixed_dt
                           << "; max. CFL from inputs = " << m_max_cfl << '\n';
        }
    }

    // Toggle initialization state and enter evolution phase
    m_is_init = false;
    if (continue_sim) {
        for (int i = max_time_states - 1; i > 0; --i) {
            m_dt[i] = m_dt[i - 1];
        }

        m_time_index++;
        m_cur_time = m_new_time;
        // New time calculated in advance_time, after dt has been computed

        // clang-format off
        if (m_verbose >= 0) {
            amrex::Print()
                << "\n==============================================================================\n";
}
        // clang-format on
    } else {
        m_cur_time = m_new_time;
    }

    return continue_sim;
}

void SimTime::set_current_cfl(
    const amrex::Real conv_cfl,
    const amrex::Real diff_cfl,
    const amrex::Real src_cfl)
{
    bool use_init_dt{false};
    const amrex::Real cd_cfl = conv_cfl + diff_cfl;
    const amrex::Real cfl_unit_time =
        cd_cfl + std::sqrt((cd_cfl * cd_cfl) + (4.0_rt * src_cfl));
    if ((m_adaptive && !m_is_init) &&
        (cfl_unit_time < std::numeric_limits<amrex::Real>::epsilon())) {
        // First timestep, starting from t = 0, is special case
        if (m_cur_time < std::numeric_limits<amrex::Real>::epsilon()) {
            if (m_initial_dt > 0.0_rt) {
                use_init_dt = true;
            } else {
                amrex::Abort(
                    "CFL is below machine epsilon, the time step is adaptive, "
                    "and the initial dt is not set. Please set the initial dt, "
                    "use a fixed time step, or fix the case setup");
            }
        } else {
            amrex::Abort(
                "CFL is below machine epsilon and the time step is adaptive. "
                "Please use a fixed time step or fix the case setup");
        }
    }
    amrex::Real dt_new =
        use_init_dt ? m_initial_dt
                    : 2.0_rt * m_max_cfl /
                          amrex::max(
                              cfl_unit_time,
                              std::numeric_limits<amrex::Real>::epsilon());

    // Restrict timestep during initialization phase
    if (m_is_init) {
        dt_new *= m_init_shrink;
    }

    // Limit timestep growth to % per timestep
    if (m_dt[0] > 0.0_rt && m_adaptive) {
        dt_new =
            amrex::min<amrex::Real>(dt_new, (1.0_rt + m_dt_growth) * m_dt[0]);
    }

    if (m_adaptive) {
        if (m_max_dt > 0.0_rt) {
            dt_new = amrex::min(dt_new, m_max_dt);
        }

        if (m_min_dt > 0.0_rt && dt_new < m_min_dt) {
            amrex::Abort(
                "CFL restriction on adaptive dt has reduced the time step size "
                "below the minimum allowable dt (min_dt). Aborting simulation "
                "because velocity values have likely become too large and "
                "unphysical.");
        }

        // Don't overshoot stop time
        if ((m_stop_time > 0.0_rt) && ((m_cur_time + dt_new) > m_stop_time)) {
            dt_new = m_stop_time - m_cur_time;
        }

        // Shorten timestep to hit output frequency exactly
        // Only should be active after delay interval is passed
        if (m_chkpt_t_interval > 0.0_rt && m_force_chkpt_dt &&
            m_cur_time + dt_new - m_chkpt_t_delay >= 0) {
            // Shorten dt if going to overshoot next output time
            dt_new = get_enforced_dt_for_output(
                dt_new, m_cur_time, m_chkpt_t_interval, m_force_chkpt_tol);
            // how it works: the floor operator gets the index of the last
            // output, with a tolerance proportional to the current dt.
            // adding 1 and multiplying by the time interval finds the next
            // expected output time. if the time increment between the current
            // time and the next expected output is less than the current dt,
            // then the dt becomes this increment to get to the next output
            // time.
        }
        if (m_plt_t_interval > 0.0_rt && m_force_plt_dt &&
            m_cur_time + dt_new - m_plt_t_delay >= 0) {
            // Shorten dt if going to overshoot next output time
            dt_new = get_enforced_dt_for_output(
                dt_new, m_cur_time, m_plt_t_interval, m_force_plt_tol);
        }
        // Consider postprocessing as well
        for (int npp = 0; npp < m_postprocess_time_interval.size(); ++npp) {
            if (m_postprocess_time_interval[npp] > 0.0_rt &&
                m_postprocess_enforce_dt[npp] &&
                m_cur_time + dt_new - m_postprocess_time_delay[npp] >= 0) {
                // Shorten dt if going to overshoot next output time
                dt_new = get_enforced_dt_for_output(
                    dt_new, m_cur_time, m_postprocess_time_interval[npp],
                    m_postprocess_enforce_dt_tol[npp]);
            }
        }

        if (m_is_init && m_initial_dt > 0.0_rt) {
            dt_new = amrex::min(dt_new, m_initial_dt);
        }

        m_dt[0] = (m_cur_time <
                   (m_delay_time - 0.1 * amrex::max(m_initial_dt, 0.0_rt)))
                      ? m_initial_dt
                      : dt_new;

    } else {
        // Ensure that we use user-specified dt. Checkpoint restart might have
        // overridden this
        m_dt[0] = (m_cur_time <
                   (m_delay_time - 0.1 * amrex::max(m_initial_dt, 0.0_rt)))
                      ? m_initial_dt
                      : m_fixed_dt;
    }

    m_current_cfl = 0.5_rt * cfl_unit_time * m_dt[0];
    m_conv_cfl = conv_cfl * m_dt[0];
    m_diff_cfl = diff_cfl * m_dt[0];
    m_src_cfl = std::sqrt(src_cfl) * m_dt[0];
    m_dt_calc = dt_new;

    // Print statements for during initialization
    if (m_is_init && m_verbose >= 0) {
        amrex::Print() << "dt: " << std::setprecision(6) << m_dt[0] << '\n';
        amrex::Print() << "CFL: " << std::setprecision(6) << m_current_cfl
                       << " (conv: " << std::setprecision(6) << m_conv_cfl
                       << " diff: " << std::setprecision(6) << m_diff_cfl
                       << " src: " << std::setprecision(6) << m_src_cfl << " )"
                       << '\n';
    }
}

void SimTime::advance_time()
{
    m_new_time = m_cur_time + m_dt[0];
    // In overset simulations, dt can change after CFL calculation
    // Use ratio between calculated or set dt and current dt
    amrex::Real factor = m_dt[0] / (m_adaptive ? m_dt_calc : m_fixed_dt);
    if (m_verbose >= 0) {
        amrex::Print() << "Step: " << m_time_index << " dt: " << m_dt[0]
                       << " Time: " << m_cur_time << " to " << m_new_time
                       << '\n';
        amrex::Print() << "CFL: " << std::setprecision(6)
                       << m_current_cfl * factor
                       << " (conv: " << std::setprecision(6)
                       << m_conv_cfl * factor
                       << " diff: " << std::setprecision(6)
                       << m_diff_cfl * factor
                       << " src: " << std::setprecision(6) << m_src_cfl * factor
                       << " )" << '\n';
    }
    // If user has specified fixed delta_t then issue a warning if the timestep
    // is larger than the delta_t determined from max. CFL considerations.
    // Only issue warnings when the error is greater than 1% of the timestep
    // specified
    const bool issue_cfl_warning =
        m_adaptive ? false : (1.0_rt - (m_dt_calc / m_dt[0])) > 0.01_rt;
    if (issue_cfl_warning && !m_is_init) {
        amrex::Print() << "WARNING: fixed_dt does not satisfy CFL condition.\n"
                       << "Max. CFL: " << m_max_cfl
                       << " => dt: " << std::setprecision(6) << m_dt_calc
                       << "; dt_inp: " << m_fixed_dt << '\n';
    }
}

bool SimTime::continue_simulation() const
{
    constexpr amrex::Real eps =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    bool stop_simulation = false;

    if (m_stop_time_index == 0) {
        return stop_simulation;
    }

    if ((m_stop_time > 0.0_rt) && ((m_new_time + eps) >= m_stop_time)) {
        return stop_simulation;
    }

    if ((m_stop_time_index > 0) && (m_time_index >= m_stop_time_index)) {
        return stop_simulation;
    }

    return !(stop_simulation);
}

bool SimTime::do_regrid() const
{
    return (
        (m_regrid_interval > 0) &&
        ((m_time_index - m_regrid_start_index) > 0) &&
        ((m_time_index - m_regrid_start_index) % m_regrid_interval == 0));
}

bool SimTime::write_plot_file() const
{
    // If dt is enforced, allow smallest tolerance to be in effect. This avoids
    // unintentionally plotting in consecutive timesteps because of shortened dt
    const amrex::Real tol =
        amrex::min<amrex::Real>(m_plt_t_tol * m_dt[0], m_force_dt_abs_tol);
    return (
        ((m_plt_interval > 0) && (m_time_index - m_plt_delay >= 0) &&
         ((m_time_index - m_plt_start_index) % m_plt_interval == 0)) ||
        ((m_plt_t_interval > 0.0_rt) &&
         (m_new_time + tol - m_plt_t_delay >= 0.0_rt) &&
         (((m_new_time - m_plt_start_time + tol) / m_plt_t_interval) -
              std::floor(
                  (m_new_time - m_plt_start_time + tol) / m_plt_t_interval) <
          m_dt[0] / m_plt_t_interval)));
}

bool SimTime::write_checkpoint() const
{
    // If dt is enforced, use smallest tolerance
    const amrex::Real tol =
        amrex::min<amrex::Real>(m_chkpt_t_tol * m_dt[0], m_force_dt_abs_tol);
    return (
        ((m_chkpt_interval > 0) && (m_time_index - m_chkpt_delay >= 0) &&
         ((m_time_index - m_chkpt_start_index) % m_chkpt_interval == 0)) ||
        ((m_chkpt_t_interval > 0.0_rt) &&
         (m_new_time + tol - m_chkpt_t_delay >= 0.0_rt) &&
         (((m_new_time - m_chkpt_start_time + tol) / m_chkpt_t_interval) -
              std::floor(
                  (m_new_time - m_chkpt_start_time + tol) /
                  m_chkpt_t_interval) <
          m_dt[0] / m_chkpt_t_interval)));
}

bool SimTime::write_last_plot_file() const
{
    const amrex::Real tol =
        amrex::min<amrex::Real>(m_plt_t_tol * m_dt[0], m_force_dt_abs_tol);
    return (
        ((m_plt_interval > 0) &&
         ((m_time_index - m_plt_start_index) % m_plt_interval != 0)) ||
        ((m_plt_t_interval > 0.0_rt) &&
         (((m_new_time - m_plt_start_time + tol) / m_plt_t_interval) -
              std::floor(
                  (m_new_time - m_plt_start_time + tol) / m_plt_t_interval) >=
          m_dt[0] / m_plt_t_interval)));
}

bool SimTime::write_last_checkpoint() const
{
    const amrex::Real tol =
        amrex::min<amrex::Real>(m_chkpt_t_tol * m_dt[0], m_force_dt_abs_tol);
    return (
        ((m_chkpt_interval > 0) &&
         ((m_time_index - m_chkpt_start_index) % m_chkpt_interval != 0)) ||
        ((m_chkpt_t_interval > 0.0_rt) &&
         (((m_new_time - m_chkpt_start_time + tol) / m_chkpt_t_interval) -
              std::floor(
                  (m_new_time - m_chkpt_start_time + tol) /
                  m_chkpt_t_interval) >=
          m_dt[0] / m_chkpt_t_interval)));
}

void SimTime::set_restart_time(int tidx, amrex::Real time)
{
    m_time_index = tidx;
    m_start_time_index = tidx;
    m_regrid_start_index = tidx;
    m_plt_start_index = tidx;
    m_chkpt_start_index = tidx;
    m_plt_start_time = time;
    m_chkpt_start_time = time;

    // check if the user has set plt/regrid/checkpoint index
    amrex::ParmParse pp("time");
    pp.query("regrid_start", m_regrid_start_index);
    pp.query("plot_start", m_plt_start_index);
    pp.query("checkpoint_start", m_chkpt_start_index);
    pp.query("plot_start_time", m_plt_start_time);
    pp.query("checkpoint_start_time", m_chkpt_start_time);

    m_new_time = time;
    m_cur_time = time;
    m_start_time = time;
}

void SimTime::calculate_minimum_enforce_dt_abs_tol()
{
    m_force_dt_abs_tol =
        (m_force_chkpt_dt
             ? amrex::min<amrex::Real>(
                   m_force_dt_abs_tol, m_force_chkpt_tol * m_chkpt_t_interval)
             : m_force_dt_abs_tol);
    m_force_dt_abs_tol =
        (m_force_plt_dt
             ? amrex::min<amrex::Real>(
                   m_force_dt_abs_tol, m_force_plt_tol * m_plt_t_interval)
             : m_force_dt_abs_tol);
    for (int npp = 0; npp < m_postprocess_time_interval.size(); ++npp) {
        m_force_dt_abs_tol =
            (m_postprocess_enforce_dt[npp]
                 ? amrex::min<amrex::Real>(
                       m_force_dt_abs_tol, m_postprocess_enforce_dt_tol[npp] *
                                               m_postprocess_time_interval[npp])
                 : m_force_dt_abs_tol);
    }
}

bool SimTime::output_profiling_info() const
{
    return (
        (m_profiling_interval > 0) &&
        (m_time_index % m_profiling_interval == 0));
}

} // namespace amr_wind
