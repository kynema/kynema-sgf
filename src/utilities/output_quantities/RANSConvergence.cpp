#include "src/utilities/output_quantities/RANSConvergence.H"
#include "src/utilities/IOManager.H"
#include "src/utilities/constants.H"
#include "src/turbulence/TurbulenceModel.H"

#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_REAL.H"

#include <cmath>
#include <fstream>
#include <limits>
#include <iomanip>
#include <utility>

using namespace amrex::literals;

namespace kynema_sgf::rans_convergence {

namespace {

//! Offsets of the monitored quantities within one window entry
constexpr int nvals_per_point = 3;
constexpr int ival_u = 0;
constexpr int ival_v = 1;
constexpr int ival_tke = 2;

} // namespace

RANSConvergence::RANSConvergence(CFDSim& sim, std::string label)
    : m_sim(sim), m_label(std::move(label))
{}

RANSConvergence::~RANSConvergence() = default;

void RANSConvergence::check_turbulence_model() const
{
    const std::string model = m_sim.turbulence_model().model_name();
    if (model != "KLAxell") {
        amrex::Abort(
            "RANSConvergence: this convergence monitor is only valid for the "
            "KLAxell RANS model, but the active turbulence model is '" +
            model +
            "'. The criterion is a drift test on an instantaneous value, "
            "which is only meaningful when that value approaches a limit. A "
            "large eddy simulation has no such limit, so the monitor would "
            "report a convergence that the solution does not have.");
    }
}

void RANSConvergence::initialize()
{
    BL_PROFILE("kynema-sgf::RANSConvergence::initialize");

    check_turbulence_model();

    {
        amrex::ParmParse pp(m_label);
        populate_output_parameters(pp);

        pp.get("start_time", m_start_time);
        pp.get("sample_interval_time", m_sample_interval);
        pp.get("window", m_window);

        pp.query("velocity_abs_tol", m_vel_abs_tol);
        pp.query("velocity_rel_tol", m_vel_rel_tol);
        pp.query("tke_abs_tol", m_tke_abs_tol);
        pp.query("tke_rel_tol", m_tke_rel_tol);
        pp.query("min_samples", m_min_samples);
        pp.query("hold_time", m_hold_time);
        pp.query("stop_on_convergence", m_stop_on_convergence);

        pp.query("report_eta", m_report_eta);
        pp.query("eta_fit_window", m_eta_fit_window);
        pp.query("eta_min_samples", m_eta_min_samples);
        pp.query("eta_min_rsq", m_eta_min_rsq);
    }

    if (m_start_time <= 0.0_rt) {
        amrex::Abort(
            "RANSConvergence: " + m_label +
            ".start_time must be greater than zero. A KLAxell run is driven "
            "toward its target by forcing terms that ramp in over time, so "
            "convergence measured before those ramps have finished says only "
            "that the solution is tracking a moving target.");
    }
    if (m_sample_interval <= 0.0_rt) {
        amrex::Abort(
            "RANSConvergence: " + m_label +
            ".sample_interval_time must be greater than zero.");
    }
    if (m_window <= m_sample_interval) {
        amrex::Abort(
            "RANSConvergence: " + m_label +
            ".window must be larger than the sample interval, otherwise the "
            "window holds a single sample and its spread is always zero.");
    }
    // Every convergence test divides the spread by this tolerance, so it has
    // to be strictly positive. A zero absolute tolerance would divide by zero
    // once the relative term vanished, which raises under the FPE traps the
    // tests run with, and it is unsatisfiable in any case: the spread would
    // have to be exactly zero
    if (m_vel_abs_tol <= 0.0_rt || m_tke_abs_tol <= 0.0_rt) {
        amrex::Abort(
            "RANSConvergence: " + m_label + ".velocity_abs_tol and " + m_label +
            ".tke_abs_tol must both be greater than zero. They are the floor "
            "that keeps the test meaningful where the monitored quantity is "
            "small; a value of zero can never be met.");
    }
    if (m_vel_rel_tol < 0.0_rt || m_tke_rel_tol < 0.0_rt) {
        amrex::Abort(
            "RANSConvergence: " + m_label + ".velocity_rel_tol and " + m_label +
            ".tke_rel_tol must not be negative.");
    }
    if (m_min_samples < 2) {
        amrex::Abort(
            "RANSConvergence: " + m_label +
            ".min_samples must be at least 2 for a spread to be meaningful.");
    }
    if (m_hold_time < 0.0_rt) {
        // Two windows. The envelope has its own dip at each turning point of
        // the inertial oscillation, of depth set by the ratio of the window to
        // the oscillation period, and that dip lasts on the order of one
        // window. Requiring the criterion to hold for longer than the dip
        // means a turning point cannot stop the run however tight the window.
        m_hold_time = 2.0_rt * m_window;
    }
    if (m_eta_fit_window <= 0.0_rt) {
        // Default to a fit history several windows long, so that the fit sees
        // an actual trend rather than the noise within one window
        m_eta_fit_window = 10.0_rt * m_window;
    }

    setup_container();
    reject_terrain_points();
    prepare_ascii_file();

    m_next_sample_time =
        amrex::max<amrex::Real>(m_start_time, m_sim.time().new_time());

    amrex::Print() << "RANSConvergence: monitoring " << m_npts
                   << " points, starting at t = " << m_start_time
                   << " s, sampling every " << m_sample_interval << " s over a "
                   << m_window << " s window\n";
    for (int i = 0; i < m_npts; ++i) {
        amrex::Print() << "  point " << i << ": ("
                       << m_point_coords[(AMREX_SPACEDIM * i) + 0] << ", "
                       << m_point_coords[(AMREX_SPACEDIM * i) + 1] << ", "
                       << m_point_coords[(AMREX_SPACEDIM * i) + 2] << ")\n";
    }
    if (!m_stop_on_convergence) {
        amrex::Print() << "RANSConvergence: reporting only, the run will not "
                          "be stopped by this monitor\n";
    }
}

void RANSConvergence::setup_container()
{
    BL_PROFILE("kynema-sgf::RANSConvergence::setup_container");

    auto& repo = m_sim.repo();

    m_fields.emplace_back(&repo.get_field("velocity"));
    // KLAxell always solves a tke transport equation, so the field exists
    m_fields.emplace_back(&repo.get_field("tke"));

    int ncomp = 0;
    for (const auto* fld : m_fields) {
        ncomp += fld->num_comp();
    }

    m_has_terrain = repo.field_exists("terrain_height");
    m_ncomp_total = ncomp + (m_has_terrain ? 1 : 0);

    m_sampler = sampling::SamplerBase::create("ProbeSampler", m_sim);
    m_sampler->label() = m_label;
    m_sampler->id() = 0;
    m_sampler->initialize(m_label);

    m_npts = static_cast<int>(m_sampler->num_points());
    if (m_npts < 1) {
        amrex::Abort(
            "RANSConvergence: no monitor points were read for " + m_label);
    }

    sampling::SampleLocType locs;
    m_sampler->sampling_locations(locs);
    m_point_coords.resize(static_cast<size_t>(AMREX_SPACEDIM) * m_npts);
    for (int i = 0; i < m_npts; ++i) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            m_point_coords[(AMREX_SPACEDIM * i) + d] = locs.locations()[i][d];
        }
    }

    amrex::Vector<std::unique_ptr<sampling::SamplerBase>> samplers;
    samplers.emplace_back(std::move(m_sampler));

    m_scontainer = std::make_unique<sampling::SamplingContainer>(m_sim.mesh());
    m_scontainer->setup_container(m_ncomp_total);
    m_scontainer->set_interpolation_order(1);
    m_scontainer->initialize_particles(samplers);
    m_scontainer->Redistribute();
    m_scontainer->num_sampling_particles() = m_npts;

    m_sampler = std::move(samplers[0]);

    m_sample_buf.assign(static_cast<size_t>(m_ncomp_total) * m_npts, 0.0_rt);
}

void RANSConvergence::reject_terrain_points()
{
    BL_PROFILE("kynema-sgf::RANSConvergence::reject_terrain_points");

    if (!m_has_terrain) {
        return;
    }

    // Interpolate the terrain surface height once, into the spare component
    const int ncomp_flow = m_ncomp_total - 1;
    amrex::Vector<Field*> terrain_fields;
    terrain_fields.emplace_back(&m_sim.repo().get_field("terrain_height"));
    m_scontainer->interpolate_fields(terrain_fields, ncomp_flow);
    m_scontainer->populate_buffer(m_sample_buf);

    int nbad = 0;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        const long offset = static_cast<long>(ncomp_flow) * m_npts;
        for (int i = 0; i < m_npts; ++i) {
            const amrex::Real zterrain = m_sample_buf[offset + i];
            const amrex::Real zpt = m_point_coords[(AMREX_SPACEDIM * i) + 2];
            if (zpt < zterrain) {
                ++nbad;
                amrex::Print()
                    << "RANSConvergence: monitor point " << i << " at ("
                    << m_point_coords[(AMREX_SPACEDIM * i) + 0] << ", "
                    << m_point_coords[(AMREX_SPACEDIM * i) + 1] << ", " << zpt
                    << ") lies below the terrain surface at z = " << zterrain
                    << '\n';
            }
        }
    }
    amrex::ParallelDescriptor::ReduceIntMax(nbad);

    if (nbad > 0) {
        amrex::Abort(
            "RANSConvergence: monitor points lie inside the terrain. Such a "
            "point is held near zero by the drag forcing and converges "
            "immediately, and because every point must pass the test it would "
            "weaken the criterion without ever failing. Move the points listed "
            "above into the fluid.");
    }
}

void RANSConvergence::post_regrid_actions()
{
    BL_PROFILE("kynema-sgf::RANSConvergence::post_regrid_actions");
    if (m_scontainer) {
        m_scontainer->Redistribute();
    }
}

void RANSConvergence::post_advance_work()
{
    BL_PROFILE("kynema-sgf::RANSConvergence::post_advance_work");

    if (m_converged) {
        return;
    }

    const amrex::Real cur_time = m_sim.time().new_time();
    if (cur_time < m_start_time) {
        return;
    }
    // The cadence is decided from replicated state, so every rank enters the
    // collective sampling calls on the same step
    if (cur_time < m_next_sample_time) {
        return;
    }

    take_sample();
    evaluate_convergence();

    m_next_sample_time = cur_time + m_sample_interval;
}

void RANSConvergence::take_sample()
{
    BL_PROFILE("kynema-sgf::RANSConvergence::take_sample");

    m_scontainer->interpolate_fields(m_fields, 0);
    m_scontainer->populate_buffer(m_sample_buf);

    const amrex::Real cur_time = m_sim.time().new_time();

    // Only the IO rank holds a fully reduced buffer, and it is the only rank
    // that keeps the window, computes the spreads and decides
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    amrex::Vector<amrex::Real> entry(
        static_cast<size_t>(nvals_per_point) * m_npts, 0.0_rt);
    const long np = m_npts;
    for (int i = 0; i < m_npts; ++i) {
        entry[(nvals_per_point * i) + ival_u] = m_sample_buf[i];
        entry[(nvals_per_point * i) + ival_v] = m_sample_buf[np + i];
        // Component 2 is the vertical velocity, which is deliberately not
        // monitored; component 3 is tke
        entry[(nvals_per_point * i) + ival_tke] = m_sample_buf[(3 * np) + i];
    }

    m_times.push_back(cur_time);
    m_samples.push_back(std::move(entry));

    // Drop the entries that have aged out of the window. The window holds a
    // few dozen entries at most, so erasing from the front is cheap
    while (m_times.size() > 1 && (cur_time - m_times.front()) > m_window) {
        m_times.erase(m_times.begin());
        m_samples.erase(m_samples.begin());
    }
}

EnvelopeFit RANSConvergence::fit_envelope_decay(
    const amrex::Vector<amrex::Real>& times,
    const amrex::Vector<amrex::Real>& spreads,
    const amrex::Real threshold,
    const int min_samples)
{
    EnvelopeFit fit;

    if (times.size() != spreads.size()) {
        return fit;
    }
    const int n = static_cast<int>(times.size());
    if (n < min_samples || n < 2) {
        return fit;
    }
    if (threshold <= 0.0_rt) {
        return fit;
    }
    // A zero or negative spread has no logarithm. It also means the envelope
    // has already collapsed, so there is nothing left to extrapolate
    for (const auto spread : spreads) {
        if (spread <= 0.0_rt) {
            return fit;
        }
    }

    // Ordinary least squares of ln(s) against t, which is the model
    // s(t) = A exp(-rate * t)
    const auto rn = static_cast<amrex::Real>(n);
    amrex::Real sum_t = 0.0_rt;
    amrex::Real sum_y = 0.0_rt;
    for (int i = 0; i < n; ++i) {
        sum_t += times[i];
        sum_y += std::log(spreads[i]);
    }
    const amrex::Real mean_t = sum_t / rn;
    const amrex::Real mean_y = sum_y / rn;

    amrex::Real stt = 0.0_rt;
    amrex::Real sty = 0.0_rt;
    amrex::Real syy = 0.0_rt;
    for (int i = 0; i < n; ++i) {
        const amrex::Real dt = times[i] - mean_t;
        const amrex::Real dy = std::log(spreads[i]) - mean_y;
        stt += dt * dt;
        sty += dt * dy;
        syy += dy * dy;
    }
    if (stt <= 0.0_rt) {
        return fit;
    }

    const amrex::Real slope = sty / stt;
    const amrex::Real intercept = mean_y - (slope * mean_t);

    fit.rate = -slope;
    fit.amplitude = std::exp(intercept);
    fit.rsq = (syy > 0.0_rt) ? ((sty * sty) / (stt * syy)) : 0.0_rt;

    // Only a decaying envelope can reach the threshold from above
    if (fit.rate <= 0.0_rt) {
        return fit;
    }
    const amrex::Real t_cross = std::log(fit.amplitude / threshold) / fit.rate;
    const amrex::Real eta = t_cross - times[n - 1];
    if (eta <= 0.0_rt) {
        // The fit says the threshold should already have been met; the
        // measured spread disagrees, so the extrapolation is not usable
        return fit;
    }

    fit.time_to_threshold = eta;
    fit.valid = true;
    return fit;
}

void RANSConvergence::evaluate_convergence()
{
    BL_PROFILE("kynema-sgf::RANSConvergence::evaluate_convergence");

    int converged_flag = 0;

    if (amrex::ParallelDescriptor::IOProcessor()) {
        const int nsamples = static_cast<int>(m_times.size());
        const bool window_full = (nsamples >= m_min_samples) &&
                                 ((m_times.back() - m_times.front()) >=
                                  (m_window - (0.5_rt * m_sample_interval)));

        int num_converged = 0;
        int worst_vel_point = -1;
        int worst_tke_point = -1;
        amrex::Real worst_vel_ratio = -1.0_rt;
        amrex::Real worst_tke_ratio = -1.0_rt;
        amrex::Real worst_vel_spread = 0.0_rt;
        amrex::Real worst_tke_spread = 0.0_rt;
        amrex::Real worst_vel_tol = 0.0_rt;
        amrex::Real worst_tke_tol = 0.0_rt;

        for (int i = 0; i < m_npts; ++i) {
            amrex::Real smin = std::numeric_limits<amrex::Real>::max();
            amrex::Real smax = std::numeric_limits<amrex::Real>::lowest();
            amrex::Real ssum = 0.0_rt;
            amrex::Real kmin = std::numeric_limits<amrex::Real>::max();
            amrex::Real kmax = std::numeric_limits<amrex::Real>::lowest();
            amrex::Real ksum = 0.0_rt;

            for (const auto& entry : m_samples) {
                const amrex::Real uu = entry[(nvals_per_point * i) + ival_u];
                const amrex::Real vv = entry[(nvals_per_point * i) + ival_v];
                // Horizontal speed only: the vertical component is small and
                // comparatively noisy in a converging boundary layer, and
                // including it would add jitter to the envelope
                const amrex::Real spd = std::sqrt((uu * uu) + (vv * vv));
                const amrex::Real kk = entry[(nvals_per_point * i) + ival_tke];

                smin = amrex::min<amrex::Real>(smin, spd);
                smax = amrex::max<amrex::Real>(smax, spd);
                ssum += spd;
                kmin = amrex::min<amrex::Real>(kmin, kk);
                kmax = amrex::max<amrex::Real>(kmax, kk);
                ksum += kk;
            }

            const auto rn = static_cast<amrex::Real>(m_samples.size());
            const amrex::Real vel_spread = smax - smin;
            const amrex::Real tke_spread = kmax - kmin;
            const amrex::Real vel_tol =
                effective_tolerance(ssum / rn, m_vel_abs_tol, m_vel_rel_tol);
            const amrex::Real tke_tol =
                effective_tolerance(ksum / rn, m_tke_abs_tol, m_tke_rel_tol);

            const amrex::Real vel_ratio = vel_spread / vel_tol;
            const amrex::Real tke_ratio = tke_spread / tke_tol;

            if (vel_ratio > worst_vel_ratio) {
                worst_vel_ratio = vel_ratio;
                worst_vel_point = i;
                worst_vel_spread = vel_spread;
                worst_vel_tol = vel_tol;
            }
            if (tke_ratio > worst_tke_ratio) {
                worst_tke_ratio = tke_ratio;
                worst_tke_point = i;
                worst_tke_spread = tke_spread;
                worst_tke_tol = tke_tol;
            }
            if (vel_ratio < 1.0_rt && tke_ratio < 1.0_rt) {
                ++num_converged;
            }
        }

        // The normalized spreads are what must fall below one, so they are
        // the series the exponential fit extrapolates
        // Only full windows enter the fit history. While the window is
        // filling it holds a couple of samples, whose spread is small whatever
        // the flow is doing, and those points would flatten the fitted decay
        const amrex::Real cur_time = m_times.back();
        if (window_full) {
            m_eta_times.push_back(cur_time);
            m_eta_vel_ratio.push_back(worst_vel_ratio);
            m_eta_tke_ratio.push_back(worst_tke_ratio);
            while (m_eta_times.size() > 1 &&
                   (cur_time - m_eta_times.front()) > m_eta_fit_window) {
                m_eta_times.erase(m_eta_times.begin());
                m_eta_vel_ratio.erase(m_eta_vel_ratio.begin());
                m_eta_tke_ratio.erase(m_eta_tke_ratio.begin());
            }
        }

        amrex::Real eta = -1.0_rt;
        amrex::Real eta_rate = 0.0_rt;
        amrex::Real eta_rsq = 0.0_rt;
        if (m_report_eta && window_full) {
            const auto vfit = fit_envelope_decay(
                m_eta_times, m_eta_vel_ratio, 1.0_rt, m_eta_min_samples);
            const auto kfit = fit_envelope_decay(
                m_eta_times, m_eta_tke_ratio, 1.0_rt, m_eta_min_samples);

            // Both quantities must converge, so the later estimate governs
            if (vfit.valid && vfit.rsq >= m_eta_min_rsq) {
                eta = vfit.time_to_threshold;
                eta_rate = vfit.rate;
                eta_rsq = vfit.rsq;
            }
            if (kfit.valid && kfit.rsq >= m_eta_min_rsq &&
                kfit.time_to_threshold > eta) {
                eta = kfit.time_to_threshold;
                eta_rate = kfit.rate;
                eta_rsq = kfit.rsq;
            }
        }

        amrex::Print() << "RANSConvergence: t = " << std::scientific
                       << std::setprecision(4) << cur_time << " s, converged "
                       << num_converged << "/" << m_npts << " points";
        if (!window_full) {
            amrex::Print() << " (filling window, " << nsamples << " samples)";
        }
        amrex::Print() << '\n';
        if (worst_vel_point >= 0) {
            amrex::Print()
                << "  worst speed: point " << worst_vel_point << " ("
                << m_point_coords[(AMREX_SPACEDIM * worst_vel_point) + 0]
                << ", "
                << m_point_coords[(AMREX_SPACEDIM * worst_vel_point) + 1]
                << ", "
                << m_point_coords[(AMREX_SPACEDIM * worst_vel_point) + 2]
                << ") spread " << worst_vel_spread << " m/s, tol "
                << worst_vel_tol << " m/s\n";
        }
        if (worst_tke_point >= 0) {
            amrex::Print()
                << "  worst tke  : point " << worst_tke_point << " ("
                << m_point_coords[(AMREX_SPACEDIM * worst_tke_point) + 0]
                << ", "
                << m_point_coords[(AMREX_SPACEDIM * worst_tke_point) + 1]
                << ", "
                << m_point_coords[(AMREX_SPACEDIM * worst_tke_point) + 2]
                << ") spread " << worst_tke_spread << ", tol " << worst_tke_tol
                << '\n';
        }
        if (eta > 0.0_rt) {
            amrex::Print() << "  estimated time to convergence: " << eta
                           << " s (decay rate " << eta_rate
                           << " 1/s, R2 = " << std::defaultfloat
                           << std::setprecision(3) << eta_rsq << ")\n";
        }

        // Passing once is not enough. The envelope dips at every turning
        // point of the inertial oscillation, so the criterion must hold
        // continuously for longer than such a dip lasts
        const bool all_pass = window_full && (num_converged == m_npts);
        if (all_pass) {
            if (m_converged_since < 0.0_rt) {
                m_converged_since = cur_time;
            }
        } else {
            m_converged_since = -1.0_rt;
        }
        const amrex::Real hold_elapsed =
            all_pass ? (cur_time - m_converged_since) : 0.0_rt;

        if (all_pass && hold_elapsed < m_hold_time) {
            amrex::Print() << "  all points within tolerance, holding for "
                           << hold_elapsed << " of " << m_hold_time << " s\n";
        }

        write_ascii(
            nsamples, window_full, hold_elapsed, num_converged, worst_vel_point,
            worst_vel_spread, worst_vel_tol, worst_tke_point, worst_tke_spread,
            worst_tke_tol, eta);

        if (all_pass && hold_elapsed >= m_hold_time) {
            converged_flag = 1;
        }
    }

    // Every rank must reach the same verdict or they will disagree about
    // entering the next timestep and deadlock in the following collective
    amrex::ParallelDescriptor::ReduceIntMax(converged_flag);
    m_converged = (converged_flag == 1);

    if (m_converged) {
        amrex::Print() << "RANSConvergence: all " << m_npts
                       << " monitor points held within tolerance for "
                       << m_hold_time << " s\n";
        if (m_stop_on_convergence) {
            m_sim.time().request_stop(
                "RANSConvergence: all monitor points converged");
        } else {
            // Reporting only: keep watching rather than latching
            m_converged = false;
        }
    }
}

void RANSConvergence::prepare_ascii_file()
{
    BL_PROFILE("kynema-sgf::RANSConvergence::prepare_ascii_file");

    const std::string post_dir = m_sim.io_manager().post_processing_directory();
    const std::string sname =
        amrex::Concatenate(m_label, m_sim.time().time_index());
    m_out_fname = post_dir + "/" + sname + ".txt";

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream f(m_out_fname.c_str());
        // samples_in_window and window_full are recorded because
        // num_converged is meaningless while the window is still filling: a
        // window holding two samples has a small spread whatever the flow is
        // doing, so a reader of this file alone must be able to tell the two
        // apart
        f << "time samples_in_window window_full num_converged num_points "
             "worst_velocity_point worst_velocity_spread worst_velocity_tol "
             "worst_tke_point worst_tke_spread worst_tke_tol hold_elapsed "
             "estimated_time_to_convergence\n";
        f.close();
    }
}

void RANSConvergence::write_ascii(
    const int num_samples,
    const bool window_full,
    const amrex::Real hold_elapsed,
    const int num_converged,
    const int worst_vel_point,
    const amrex::Real worst_vel_spread,
    const amrex::Real worst_vel_tol,
    const int worst_tke_point,
    const amrex::Real worst_tke_spread,
    const amrex::Real worst_tke_tol,
    const amrex::Real eta)
{
    BL_PROFILE("kynema-sgf::RANSConvergence::write_ascii");

    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    std::ofstream f(m_out_fname.c_str(), std::ios_base::app);
    f << std::scientific << std::setprecision(m_precision) << std::setw(m_width)
      << m_sim.time().new_time() << ' ' << num_samples << ' '
      << (window_full ? 1 : 0) << ' ' << num_converged << ' ' << m_npts << ' '
      << worst_vel_point << std::setw(m_width) << worst_vel_spread
      << std::setw(m_width) << worst_vel_tol << ' ' << worst_tke_point
      << std::setw(m_width) << worst_tke_spread << std::setw(m_width)
      << worst_tke_tol << std::setw(m_width) << hold_elapsed
      << std::setw(m_width) << eta << '\n';
    f.close();
}

} // namespace kynema_sgf::rans_convergence
