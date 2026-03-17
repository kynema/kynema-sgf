#ifndef POSTPROCESSING_H
#define POSTPROCESSING_H

#include <memory>
#include "amr-wind/core/Factory.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

/**
 *  \defgroup utilities Utilities
 *  Utilities
 *
 *  This group documents various utilities available in AMR-Wind for
 *  input/output, post-processing, collecting statistics etc.
 */

namespace amr_wind {

class CFDSim;

/** Abstract representation of a post-processing utility
 *  \ingroup utilities
 *
 *  Provides an interface for a runtime-selectable post-processing action (e.g.,
 *  data sampling, volume/surface integration, etc.) with the main solver.
 *
 *  All post-processing utilities must derive from this class.
 */
class PostProcessBase
    : public Factory<PostProcessBase, CFDSim&, const std::string&>
{
public:
    static std::string base_identifier() { return "PostProcessBase"; }

    ~PostProcessBase() override = default;

    virtual void pre_init_actions() = 0;

    /** Initialize the post-processing utility
     *
     *  When this method is called, the mesh and the solution fields are fully
     *  initialized
     */
    virtual void initialize() = 0;

    /** Perform actions at the end of each timestep
     */
    virtual void post_advance_work() = 0;

    /** Perform output routine
     */
    virtual void output_actions() = 0;

    //! Actions to perform post regrid
    virtual void post_regrid_actions() = 0;

    void populate_output_parameters(amrex::ParmParse& pp)
    {
        pp.query("output_interval", m_out_interval);
        pp.query("output_time_interval", m_out_time_interval);
        pp.query("output_delay", m_out_delay);
        pp.query("output_time_delay", m_out_time_delay);
        pp.query("output_from_restart", m_out_from_restart);
        pp.query("output_start", m_out_start_index);
        pp.query("output_start_time", m_out_start_time);
        pp.query("enforce_output_time_dt", m_enforce_dt);
        pp.query("enforce_output_time_dt_reltol", m_enforce_dt_tol);
        pp.query("output_after_final_step", m_do_final_output);

        // Check for presence of arguments
        const bool has_interval = pp.contains("output_interval");
        const bool has_time_int = pp.contains("output_time_interval");
        const bool has_frequency = pp.contains("output_frequency");

        if (has_frequency) {
            if (has_interval) {
                amrex::Print()
                    << "WARNING: output_interval and output_frequency "
                       "(legacy argument) are both specified; "
                       "output_frequency will be ignored.\n";
            } else {
                pp.get("output_frequency", m_out_interval);
                amrex::Print()
                    << "WARNING: legacy input argument of output_frequency has "
                       "been specified; the code will function as normal, but "
                       "the new input argument of output_interval is "
                       "preferred.\n";
            }
        }

        if (m_out_from_restart &&
            (pp.contains("output_start") || pp.contains("output_start_time"))) {
            amrex::Abort(
                "Conflict in output parameters. output_from_restart is true, "
                "designating that output timing calculations consider the "
                "restart time or time step, but the output_start or "
                "output_start_time is specified, which manually overwrites the "
                "stored restart time or time step.\n");
        }

        // Because of default, interval is valid unless turned off;
        // time_interval is invalid unless turned on
        const bool has_valid_interval =
            (has_interval || has_frequency) ? m_out_interval > 0 : true;
        const bool has_valid_interval_present =
            (has_interval || has_frequency) ? m_out_interval > 0 : false;
        const bool has_valid_time_int =
            has_time_int ? m_out_time_interval > 0. : false;

        if (!(has_valid_interval || has_valid_time_int)) {
            amrex::Abort(
                "post_processing output timing parameters are specified, but "
                "none are valid. Please specify an output_interval (or "
                "output_frequency) greater than 0 or an output_time_interval "
                "greater than 0.0_rt, or rely on the default output_interval.");
        }

        if (has_valid_interval_present && has_valid_time_int) {
            amrex::Abort(
                "output_interval (or output_frequency) and "
                "output_time_interval are both specified. timestep- and "
                "time-based output parameters should not be used together; "
                "please only specify one.");
        }

        if (!has_valid_time_int && m_enforce_dt) {
            amrex::Abort(
                "enforce_output_time_dt is true, but output_time_interval "
                "has not been specified.");
        }

        // Because of default, turn off step interval when time interval is
        // meant to replace it
        if (has_valid_time_int && !(has_interval || has_frequency)) {
            m_out_interval = -1;
        }
    }

    //! Determine if output will happen this time step
    [[nodiscard]] bool do_output_now(
        const int ntime,
        const int ntime0_in,
        const amrex::Real time,
        const amrex::Real time0_in,
        const amrex::Real dt,
        amrex::Real tol) const
    {
        // Choose smallest tolerance between incoming and local
        tol = amrex::min<amrex::Real>(tol, m_out_time_tol * dt);

        // Should the restart time always be used
        int ntime0 = m_out_from_restart ? ntime0_in : 0;
        amrex::Real time0 = m_out_from_restart ? time0_in : 0_rt;
        // Has the manual start been specified
        ntime0 = m_out_start_index >= 0 ? m_out_start_index : ntime0;
        time0 = m_out_start_time > -0.5_rt ? m_out_start_time : time0;

        const bool step_based =
            (m_out_interval > 0) && ((ntime - m_out_delay >= 0) &&
                                     ((ntime - ntime0) % m_out_interval == 0));
        const bool time_based =
            (m_out_time_interval > 0.0_rt) &&
            ((time + tol - m_out_time_delay >= 0.0_rt) &&
             (time - time0 + tol -
                  (std::floor((time - time0 + tol) / m_out_time_interval) *
                   m_out_time_interval) <
              dt));

        return (step_based || time_based);
    }

    [[nodiscard]] amrex::Real output_time_interval() const
    {
        return m_out_time_interval;
    }
    [[nodiscard]] amrex::Real output_time_delay() const
    {
        return m_out_time_delay;
    }
    [[nodiscard]] amrex::Real enforce_dt_tolerance() const
    {
        return m_enforce_dt_tol;
    }
    [[nodiscard]] bool enforce_dt() const { return m_enforce_dt; }
    [[nodiscard]] bool do_final_output() const { return m_do_final_output; }

protected:
    //! Time step interval for output
    int m_out_interval{10};
    //! Time interval for output
    amrex::Real m_out_time_interval{-1.0_rt};
    //! Time step delay for output
    int m_out_delay{0};
    //! Time delay for output
    amrex::Real m_out_time_delay{0.0_rt};
    //! Always consider restart time for calculating output
    bool m_out_from_restart{false};
    //! Time step start for calculating step interval (restarts)
    int m_out_start_index{-1};
    //! Time start for calculating time interval (restarts)
    amrex::Real m_out_start_time{-1.0_rt};
    //! Relative (to dt) tolerance for time interval output
    amrex::Real m_out_time_tol{1.0e-2_rt};
    //! Relative (to m_out_time_interval) tolerance for enforcing dt
    amrex::Real m_enforce_dt_tol{1.0e-3_rt};
    //! Flag for enforcing time step based on output time
    bool m_enforce_dt{false};
    //! Flag for output at simulation end - off by default
    bool m_do_final_output{false};
};

/** A collection of post-processing instances that are active during a
 *  simulation
 *
 *  Upon initialization, the class reads the ``incflo.post_processing``
 *  parameter to determine the utilities that must be activated and uses runtime
 *  selection to initialize them
 *
 *  \ingroup utilities
 */
class PostProcessManager
{
public:
    explicit PostProcessManager(CFDSim& sim);

    ~PostProcessManager() = default;

    void pre_init_actions();

    /** Initialize post-processing utilities
     *
     *  Reads user inputs and instantiates all the registered post-processing
     *  utilities. Executes the registered utilities once before starting the
     *  time integration
     */
    void post_init_actions();

    //! Call all registered utilities to perform actions after a timestep
    void post_advance_work();

    //! Call all registered utilities to output at final time
    void final_output();

    void post_regrid_actions();

private:
    CFDSim& m_sim;

    amrex::Vector<std::unique_ptr<PostProcessBase>> m_post;
};

} // namespace amr_wind

#endif /* POSTPROCESSING_H */
