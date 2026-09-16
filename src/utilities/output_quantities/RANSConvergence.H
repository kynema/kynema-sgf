#ifndef RANSCONVERGENCE_H
#define RANSCONVERGENCE_H

#include <memory>
#include <string>

#include "src/CFDSim.H"
#include "src/utilities/PostProcessing.H"
#include "src/utilities/sampling/SamplerBase.H"
#include "src/utilities/sampling/SamplingContainer.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"

using namespace amrex::literals;

namespace kynema_sgf::rans_convergence {

/** Result of fitting an exponential decay to an envelope history
 *  \ingroup output_quantities
 *
 *  Holds the outcome of a least-squares fit of \f$s(t) = A e^{-rt}\f$ to a
 *  sequence of envelope spreads, along with the extrapolated time at which the
 *  fit would reach a requested threshold.
 */
struct EnvelopeFit
{
    //! True if the fit succeeded and the extrapolation is usable
    bool valid{false};

    //! Decay rate \f$r\f$ in 1/s; positive means the envelope is shrinking
    amrex::Real rate{0.0_rt};

    //! Amplitude \f$A\f$ of the fitted exponential
    amrex::Real amplitude{0.0_rt};

    //! Coefficient of determination of the fit, computed in log space
    amrex::Real rsq{0.0_rt};

    //! Estimated seconds from the last sample until the threshold is reached
    amrex::Real time_to_threshold{-1.0_rt};
};

/** Pseudo-steady convergence monitor for KLAxell RANS runs
 *  \ingroup output_quantities
 *
 *  Tracks horizontal wind speed and turbulent kinetic energy at a
 *  user-supplied cloud of monitor points. A point is considered converged when
 *  the peak-to-trough spread of both quantities over a trailing time window
 *  falls below tolerance. Once every point has converged the monitor can ask
 *  the time integrator to stop, which triggers a final checkpoint and
 *  plotfile.
 *
 *  The spread over a window is used in preference to a difference between
 *  successive samples because an atmospheric boundary layer under Coriolis and
 *  geostrophic forcing approaches its steady state through a damped inertial
 *  oscillation. A difference between successive samples vanishes at every
 *  turning point of that oscillation and would report convergence long before
 *  the oscillation has decayed.
 *
 *  This utility is restricted to the KLAxell turbulence model. The criterion
 *  is a drift test on an instantaneous value, which is only meaningful when
 *  that value approaches a limit; a large eddy simulation has no such limit.
 */
class RANSConvergence : public PostProcessBase::Register<RANSConvergence>
{
public:
    static std::string identifier() { return "RANSConvergence"; }

    RANSConvergence(CFDSim& /*sim*/, std::string /*label*/);

    ~RANSConvergence() override;

    //! Perform actions before mesh is created
    void pre_init_actions() override {}

    //! Read user inputs, build the monitor points and validate them
    void initialize() override;

    //! Sample the monitor points and test for convergence
    void post_advance_work() override;

    //! Output is driven internally by the sampling cadence
    void output_actions() override {}

    //! Redistribute the monitor points after the mesh changes
    void post_regrid_actions() override;

    /** Fit an exponential decay to an envelope history and extrapolate
     *
     *  Performs an ordinary least-squares fit of \f$\ln s\f$ against \f$t\f$,
     *  which corresponds to the model \f$s(t) = A e^{-rt}\f$, and solves for
     *  the time at which the fitted curve reaches \p threshold.
     *
     *  The fit is reported as invalid when there are too few samples, when any
     *  sample is not positive, when the fitted envelope is not decaying, or
     *  when the threshold has already been reached.
     *
     *  \param times Sample times in seconds, in increasing order
     *  \param spreads Envelope spread at each sample time
     *  \param threshold Value the extrapolation is solved for
     *  \param min_samples Fewest samples that will be fitted
     *  \return Fit coefficients and the extrapolated time to the threshold
     */
    static EnvelopeFit fit_envelope_decay(
        const amrex::Vector<amrex::Real>& times,
        const amrex::Vector<amrex::Real>& spreads,
        amrex::Real threshold,
        int min_samples);

    /** Tolerance actually applied to an envelope spread
     *
     *  Combines an absolute floor with a relative term so that the test stays
     *  meaningful both where the monitored quantity is small and where it is
     *  large.
     *
     *  \param mean Mean of the monitored quantity over the window
     *  \param abs_tol Absolute tolerance floor
     *  \param rel_tol Relative tolerance applied to \p mean
     *  \return The larger of the absolute floor and the relative tolerance
     */
    static amrex::Real effective_tolerance(
        const amrex::Real mean,
        const amrex::Real abs_tol,
        const amrex::Real rel_tol)
    {
        return amrex::max<amrex::Real>(abs_tol, rel_tol * std::abs(mean));
    }

private:
    //! Abort unless the active turbulence model is KLAxell
    void check_turbulence_model() const;

    //! Build the probe sampler and the particle container
    void setup_container();

    //! Abort if any monitor point lies inside the immersed terrain
    void reject_terrain_points();

    //! Interpolate the monitored fields and append them to the window
    void take_sample();

    //! Test every point against tolerance and report
    void evaluate_convergence();

    //! Prepare the ASCII output file and write its header
    void prepare_ascii_file();

    //! Append one line of convergence diagnostics
    void write_ascii(
        int num_samples,
        bool window_full,
        amrex::Real hold_elapsed,
        int num_converged,
        int worst_vel_point,
        amrex::Real worst_vel_spread,
        amrex::Real worst_vel_tol,
        int worst_tke_point,
        amrex::Real worst_tke_spread,
        amrex::Real worst_tke_tol,
        amrex::Real eta);

    //! Reference to the CFD simulation
    CFDSim& m_sim;

    //! Name of this post-processing object, used to read its inputs
    const std::string m_label;

    //! Probe sampler holding the monitor point locations
    std::unique_ptr<sampling::SamplerBase> m_sampler;

    //! Particle container performing the interpolation
    std::unique_ptr<sampling::SamplingContainer> m_scontainer;

    //! Fields interpolated at the monitor points: velocity and tke
    amrex::Vector<Field*> m_fields;

    //! Coordinates of the monitor points, gathered once at initialization
    amrex::Vector<amrex::Real> m_point_coords;

    //! Number of monitor points
    int m_npts{0};

    //! Number of real components carried by the particle container
    int m_ncomp_total{0};

    //! True when a terrain height field is available for point validation
    bool m_has_terrain{false};

    //! Buffer receiving the interpolated values, valid on the IO rank only
    amrex::Vector<amrex::Real> m_sample_buf;

    //! Simulation time at which monitoring begins
    amrex::Real m_start_time{-1.0_rt};

    //! Physical time between samples
    amrex::Real m_sample_interval{-1.0_rt};

    //! Width of the trailing window over which the envelope is measured
    amrex::Real m_window{-1.0_rt};

    //! Absolute tolerance on the horizontal wind speed envelope
    amrex::Real m_vel_abs_tol{0.01_rt};

    //! Relative tolerance on the horizontal wind speed envelope
    amrex::Real m_vel_rel_tol{0.001_rt};

    //! Absolute tolerance on the turbulent kinetic energy envelope
    amrex::Real m_tke_abs_tol{0.001_rt};

    //! Relative tolerance on the turbulent kinetic energy envelope
    amrex::Real m_tke_rel_tol{0.01_rt};

    //! Fewest samples in the window before convergence may be declared
    int m_min_samples{4};

    //! How long every point must stay within tolerance before the run stops
    amrex::Real m_hold_time{-1.0_rt};

    //! Time at which all points most recently began passing, or -1 if not
    amrex::Real m_converged_since{-1.0_rt};

    //! Ask the time integrator to stop once every point has converged
    bool m_stop_on_convergence{true};

    //! Report an estimated time to convergence from an exponential fit
    bool m_report_eta{true};

    //! Trailing span of envelope history used for the exponential fit
    amrex::Real m_eta_fit_window{-1.0_rt};

    //! Fewest envelope samples that will be fitted
    int m_eta_min_samples{5};

    //! Smallest coefficient of determination for a reported estimate
    amrex::Real m_eta_min_rsq{0.5_rt};

    //! Next simulation time at which a sample will be taken
    amrex::Real m_next_sample_time{0.0_rt};

    //! Sample times currently held in the window
    amrex::Vector<amrex::Real> m_times;

    //! Sampled values, three per point per entry: u, v and tke
    amrex::Vector<amrex::Vector<amrex::Real>> m_samples;

    //! Times at which envelope spreads were recorded, for the exponential fit
    amrex::Vector<amrex::Real> m_eta_times;

    //! History of the largest wind speed spread over all points, normalized
    //! by each point's own tolerance so that convergence is a spread of one
    amrex::Vector<amrex::Real> m_eta_vel_ratio;

    //! History of the largest tke spread over all points, normalized by each
    //! point's own tolerance so that convergence is a spread of one
    amrex::Vector<amrex::Real> m_eta_tke_ratio;

    //! True once every point has satisfied both tolerances
    bool m_converged{false};

    //! Name of the ASCII diagnostics file
    std::string m_out_fname;

    //! Field width in the ASCII output
    int m_width{18};

    //! Precision in the ASCII output
    int m_precision{10};
};

} // namespace kynema_sgf::rans_convergence

#endif /* RANSCONVERGENCE_H */
