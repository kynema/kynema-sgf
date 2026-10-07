#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>

#include "ks_test_utils/AmrexTest.H"
#include "ks_test_utils/MeshTest.H"
#include "src/CFDSim.H"
#include "src/utilities/PostProcessing.H"
#include "src/utilities/output_quantities/RANSConvergence.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

namespace rc = kynema_sgf::rans_convergence;

//! Build an exponentially decaying envelope history
amrex::Vector<amrex::Real> decaying_series(
    const amrex::Vector<amrex::Real>& times,
    const amrex::Real amplitude,
    const amrex::Real rate)
{
    amrex::Vector<amrex::Real> vals;
    vals.reserve(times.size());
    for (const auto t : times) {
        vals.push_back(amplitude * std::exp(-rate * t));
    }
    return vals;
}

amrex::Vector<amrex::Real>
uniform_times(const int n, const amrex::Real t0, const amrex::Real dt)
{
    amrex::Vector<amrex::Real> times;
    times.reserve(n);
    for (int i = 0; i < n; ++i) {
        times.push_back(t0 + (i * dt));
    }
    return times;
}

//! Absolute tolerance for comparing a value against an exact expectation.
//!
//! Scaled to the working precision so that these tests mean the same thing in
//! single and double precision. The fit takes logarithms and an exponential,
//! so a few thousand epsilon of relative error is expected rather than a few.
amrex::Real close(const amrex::Real expected)
{
    constexpr amrex::Real factor = 1.0e4_rt;
    return factor * std::numeric_limits<amrex::Real>::epsilon() *
           std::max(std::abs(expected), 1.0_rt);
}

//! Tolerance for a value read back from the monitor's ASCII diagnostics,
//! which are written with ten significant digits
amrex::Real file_close(const amrex::Real expected)
{
    return std::max(
        close(expected), 1.0e-9_rt * std::max(std::abs(expected), 1.0_rt));
}

} // namespace

TEST(RANSConvergence, effective_tolerance_takes_the_larger_term)
{
    // Relative term dominates for a large mean
    EXPECT_NEAR(
        rc::RANSConvergence::effective_tolerance(10.0_rt, 0.01_rt, 0.01_rt),
        0.1_rt, close(0.1_rt));
    // Absolute floor takes over as the mean approaches zero
    EXPECT_NEAR(
        rc::RANSConvergence::effective_tolerance(0.0_rt, 0.01_rt, 0.01_rt),
        0.01_rt, close(0.01_rt));
    // A negative mean is treated by magnitude
    EXPECT_NEAR(
        rc::RANSConvergence::effective_tolerance(-10.0_rt, 0.01_rt, 0.01_rt),
        0.1_rt, close(0.1_rt));
}

TEST(RANSConvergence, effective_tolerance_stays_positive)
{
    // Every convergence test divides a spread by this value, so it must never
    // return zero for a configuration the monitor accepts. The monitor rejects
    // a non-positive absolute tolerance at startup, which is what makes this
    // hold even where the relative term vanishes.
    EXPECT_GT(
        rc::RANSConvergence::effective_tolerance(0.0_rt, 1.0e-12_rt, 0.0_rt),
        0.0_rt);
    EXPECT_GT(
        rc::RANSConvergence::effective_tolerance(1.0e30_rt, 1.0e-12_rt, 0.0_rt),
        0.0_rt);
    // A negative relative term cannot drag the result below the floor
    EXPECT_NEAR(
        rc::RANSConvergence::effective_tolerance(10.0_rt, 0.01_rt, -1.0_rt),
        0.01_rt, close(0.01_rt));
}

TEST(RANSConvergence, envelope_fit_recovers_a_known_decay)
{
    const amrex::Real amplitude = 5.0_rt;
    const amrex::Real rate = 1.0e-3_rt;
    // Kept short enough that the series is still above the threshold at the
    // last sample, which is the only regime in which extrapolating forward
    // means anything
    const auto times = uniform_times(40, 0.0_rt, 10.0_rt);
    const auto vals = decaying_series(times, amplitude, rate);

    const auto fit =
        rc::RANSConvergence::fit_envelope_decay(times, vals, 1.0_rt, 5);

    ASSERT_TRUE(fit.valid);
    EXPECT_NEAR(fit.rate, rate, close(rate));
    EXPECT_NEAR(fit.amplitude, amplitude, close(amplitude));
    EXPECT_NEAR(fit.rsq, 1.0_rt, close(1.0_rt));

    // A exp(-rate t) = 1 gives t = ln(A)/rate, measured from the last sample
    const amrex::Real expected = (std::log(amplitude) / rate) - times.back();
    EXPECT_NEAR(fit.time_to_threshold, expected, close(expected));
}

TEST(RANSConvergence, envelope_fit_rejects_a_growing_envelope)
{
    const auto times = uniform_times(20, 0.0_rt, 100.0_rt);
    // Negative rate means the spread is growing, so it never reaches the
    // threshold from above
    const auto vals = decaying_series(times, 0.5_rt, -1.0e-3_rt);

    const auto fit =
        rc::RANSConvergence::fit_envelope_decay(times, vals, 1.0_rt, 5);

    EXPECT_FALSE(fit.valid);
    EXPECT_LT(fit.rate, 0.0_rt);
}

TEST(RANSConvergence, envelope_fit_rejects_too_few_samples)
{
    const auto times = uniform_times(4, 0.0_rt, 100.0_rt);
    const auto vals = decaying_series(times, 5.0_rt, 1.0e-3_rt);

    const auto fit =
        rc::RANSConvergence::fit_envelope_decay(times, vals, 1.0_rt, 5);

    EXPECT_FALSE(fit.valid);
}

TEST(RANSConvergence, envelope_fit_rejects_non_positive_spreads)
{
    const auto times = uniform_times(10, 0.0_rt, 100.0_rt);
    auto vals = decaying_series(times, 5.0_rt, 1.0e-3_rt);
    // An identical pair of samples collapses the spread to exactly zero,
    // which has no logarithm
    vals[4] = 0.0_rt;

    const auto fit =
        rc::RANSConvergence::fit_envelope_decay(times, vals, 1.0_rt, 5);

    EXPECT_FALSE(fit.valid);
}

TEST(RANSConvergence, envelope_fit_rejects_a_threshold_already_met)
{
    const auto times = uniform_times(20, 0.0_rt, 100.0_rt);
    // The series decays from 0.5, so it is already below a threshold of one
    const auto vals = decaying_series(times, 0.5_rt, 1.0e-3_rt);

    const auto fit =
        rc::RANSConvergence::fit_envelope_decay(times, vals, 1.0_rt, 5);

    EXPECT_FALSE(fit.valid);
    EXPECT_LT(fit.time_to_threshold, 0.0_rt);
}

TEST(RANSConvergence, envelope_fit_reports_a_poor_fit_through_rsq)
{
    // A spread that has stalled on a noise floor is not an exponential decay.
    // The fit still returns, but with an r-squared low enough that the caller
    // will discard the estimate rather than quote a confident wrong number.
    const auto times = uniform_times(20, 0.0_rt, 100.0_rt);
    amrex::Vector<amrex::Real> vals(times.size(), 2.0_rt);
    for (amrex::Long i = 0; i < vals.size(); ++i) {
        vals[i] += ((i % 2 == 0) ? 0.1_rt : -0.1_rt);
    }

    const auto fit =
        rc::RANSConvergence::fit_envelope_decay(times, vals, 1.0_rt, 5);

    EXPECT_LT(fit.rsq, 0.5_rt);
}

TEST(RANSConvergence, envelope_fit_tolerates_noise_on_a_real_decay)
{
    // The estimate does not need to be precise, only usable. A ten percent
    // multiplicative wobble should still recover the decay rate closely
    // enough to be worth printing.
    const amrex::Real amplitude = 8.0_rt;
    const amrex::Real rate = 5.0e-4_rt;
    const auto times = uniform_times(60, 0.0_rt, 20.0_rt);
    auto vals = decaying_series(times, amplitude, rate);
    for (amrex::Long i = 0; i < vals.size(); ++i) {
        vals[i] *= 1.0_rt + ((i % 3 == 0) ? 0.1_rt : -0.05_rt);
    }

    const auto fit =
        rc::RANSConvergence::fit_envelope_decay(times, vals, 1.0_rt, 5);

    ASSERT_TRUE(fit.valid);
    EXPECT_NEAR(fit.rate, rate, 0.1_rt * rate);
    // Noise of this size costs some of the fit quality, but the result stays
    // well clear of the default eta_min_rsq gate of 0.5, so the estimate
    // would be reported rather than discarded
    EXPECT_GT(fit.rsq, 0.8_rt);
}

//! Mesh-level fixture: a KLAxell setup that the monitor accepts, with the
//! velocity and tke fields set directly to scripted values each step so the
//! monitor sees a history that is known exactly
class RANSConvergenceMeshTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{8, 8, 16}};
            pp.addarr("n_cell", ncell);
            pp.add("blocking_factor", 2);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{
                {1024.0_rt, 1024.0_rt, 1024.0_rt}};
            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", 0.5_rt);
            pp.add("stop_time", 100.0_rt);
            pp.add("max_step", 14);
            pp.add("regrid_interval", -1);
            pp.add("plot_interval", -1);
            pp.add("checkpoint_interval", -1);
        }
        {
            amrex::ParmParse pp("turbulence");
            pp.add("model", m_turbulence_model);
        }
        {
            amrex::ParmParse pp("incflo");
            amrex::Vector<std::string> physics{"ABL"};
            pp.addarr("physics", physics);
            pp.add("density", 1.2_rt);
            amrex::Vector<amrex::Real> vvec{8.0_rt, 0.0_rt, 0.0_rt};
            pp.addarr("velocity", vvec);
            amrex::Vector<amrex::Real> gvec{0.0_rt, 0.0_rt, -9.81_rt};
            pp.addarr("gravity", gvec);
            pp.add("post_processing", (std::string) "convergence");
        }
        {
            amrex::ParmParse pp("ABL");
            pp.add("surface_temp_rate", 0.0_rt);
            pp.add("initial_wind_profile", true);
            amrex::Vector<amrex::Real> hts{0.0_rt, 100.0_rt, 4000.0_rt};
            pp.addarr("temperature_heights", hts);
            pp.addarr("wind_heights", hts);
            amrex::Vector<amrex::Real> t_vals{265.0_rt, 265.0_rt, 265.0_rt};
            pp.addarr("temperature_values", t_vals);
            amrex::Vector<amrex::Real> u_vals{8.0_rt, 8.0_rt, 8.0_rt};
            pp.addarr("u_values", u_vals);
            amrex::Vector<amrex::Real> v_vals{0.0_rt, 0.0_rt, 0.0_rt};
            pp.addarr("v_values", v_vals);
            amrex::Vector<amrex::Real> tke_vals{0.1_rt, 0.1_rt, 0.1_rt};
            pp.addarr("tke_values", tke_vals);
            pp.add("surface_temp_flux", 0.0_rt);
        }
        {
            amrex::ParmParse pp("transport");
            pp.add("reference_temperature", 265.0_rt);
        }
        {
            amrex::ParmParse pp("convergence");
            pp.add("type", (std::string) "RANSConvergence");
            pp.add("probe_location_file", m_probe_file);
            pp.add("start_time", 1.0_rt);
            pp.add("sample_interval_time", 0.5_rt);
            pp.add("window", m_window);
            // Low enough that the window span, not the sample count, decides
            // when the window is full
            pp.add("min_samples", m_min_samples);
            pp.add("hold_time", 1.0_rt);
            pp.add("velocity_abs_tol", 0.1_rt);
            pp.add("velocity_rel_tol", 0.0_rt);
            pp.add("tke_abs_tol", m_tke_abs_tol);
            pp.add("tke_rel_tol", 0.0_rt);
            pp.add("stop_on_convergence", m_stop_on_convergence);
            pp.add("report_eta", false);
        }
    }

    //! Two points in the fluid, at different heights
    void write_probe_file(const amrex::Real z0, const amrex::Real z1) const
    {
        std::ofstream f(m_probe_file);
        f << "2\n"
          << "512.0 512.0 " << z0 << "\n"
          << "256.0 256.0 " << z1 << "\n";
    }

    //! Build the mesh, the ABL physics and the turbulence model
    void setup_sim()
    {
        populate_parameters();
        initialize_mesh();
        sim().pde_manager().register_icns();
        sim().init_physics();
        sim().create_turbulence_model();
    }

    //! Set every cell, ghost cells included, to one value per component
    void set_fields(
        const amrex::Real u,
        const amrex::Real v,
        const amrex::Real w,
        const amrex::Real k)
    {
        auto& vel = sim().repo().get_field("velocity");
        auto& tke = sim().repo().get_field("tke");
        const int ngv = vel.num_grow()[0];
        const int ngk = tke.num_grow()[0];
        vel.setVal(u, 0, 1, ngv);
        vel.setVal(v, 1, 1, ngv);
        vel.setVal(w, 2, 1, ngv);
        tke.setVal(k, 0, 1, ngk);
    }

    void TearDown() override
    {
        remove(m_probe_file.c_str());
        remove("post_processing/convergence00000.txt");
        MeshTest::TearDown();
    }

    std::string m_turbulence_model{"KLAxell"};
    std::string m_probe_file{"rans_convergence_probes.txt"};
    amrex::Real m_window{1.5_rt};
    amrex::Real m_tke_abs_tol{0.01_rt};
    bool m_stop_on_convergence{true};
    int m_min_samples{2};
};

namespace {

//! One row of the monitor's ASCII diagnostics
struct DiagRow
{
    amrex::Real time{0.0_rt};
    int samples{0};
    int window_full{0};
    int num_converged{0};
    int num_points{0};
    int worst_vel_point{0};
    amrex::Real worst_vel_spread{0.0_rt};
    amrex::Real worst_vel_tol{0.0_rt};
    int worst_tke_point{0};
    amrex::Real worst_tke_spread{0.0_rt};
    amrex::Real worst_tke_tol{0.0_rt};
    amrex::Real hold_elapsed{0.0_rt};
    amrex::Real eta{0.0_rt};
};

amrex::Vector<DiagRow> read_diagnostics()
{
    amrex::Vector<DiagRow> rows;
    std::ifstream f("post_processing/convergence00000.txt");
    EXPECT_TRUE(f.good());
    f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    DiagRow r;
    while (f >> r.time >> r.samples >> r.window_full >> r.num_converged >>
           r.num_points >> r.worst_vel_point >> r.worst_vel_spread >>
           r.worst_vel_tol >> r.worst_tke_point >> r.worst_tke_spread >>
           r.worst_tke_tol >> r.hold_elapsed >> r.eta) {
        rows.push_back(r);
    }
    return rows;
}

//! Scripted field values at time step n (t = 0.5 n)
//!
//! The horizontal speed holds one value through step 5, drops at step 6 and
//! holds the new value from then on. The first plateau fills the window and
//! passes the test, so the hold starts; the drop must then reset it, and only
//! the second plateau may complete it. v and tke wobble within their
//! tolerances, and w swings far outside every tolerance: if the monitor ever
//! read w in place of a horizontal component or of tke, or included it in the
//! speed, no point could converge. The tke wobble is outside the tke tolerance
//! if it were read from v, so a mix-up between the two is caught as well
struct StepValues
{
    amrex::Real u, v, w, k;
};

StepValues scripted_values(const int n)
{
    const amrex::Real sgn = (n % 2 == 0) ? 1.0_rt : -1.0_rt;
    const amrex::Real u = (n <= 5) ? 12.0_rt : 10.0_rt;
    return {
        .u = u,
        .v = 3.0_rt + (0.02_rt * sgn),
        .w = 50.0_rt * sgn,
        .k = 0.5_rt + (0.001_rt * sgn)};
}

amrex::Real horizontal_speed(const StepValues& s)
{
    return std::sqrt((s.u * s.u) + (s.v * s.v));
}

//! Peak-to-trough horizontal speed spread over steps [n0, n1]
amrex::Real scripted_speed_spread(const int n0, const int n1)
{
    amrex::Real smin = std::numeric_limits<amrex::Real>::max();
    amrex::Real smax = std::numeric_limits<amrex::Real>::lowest();
    for (int n = n0; n <= n1; ++n) {
        const amrex::Real spd = horizontal_speed(scripted_values(n));
        smin = std::min(smin, spd);
        smax = std::max(smax, spd);
    }
    return smax - smin;
}

} // namespace

TEST_F(RANSConvergenceMeshTest, stops_once_the_scripted_history_settles)
{
    write_probe_file(100.0_rt, 300.0_rt);
    setup_sim();

    auto& post_manager = sim().post_manager();
    auto& time = sim().time();
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.5_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        const auto sv = scripted_values(time.time_index());
        set_fields(sv.u, sv.v, sv.w, sv.k);
        post_manager.post_advance_work();
    }

    // Sampling starts at t = 1 (step 2). The window holds four samples and is
    // full at step 5, where the first plateau passes and the hold begins. The
    // drop at step 6 resets it. The second plateau fills the window at step 9
    // (t = 4.5), and the hold is satisfied at step 11 (t = 5.5), three steps
    // before max_step
    EXPECT_TRUE(time.stop_requested());
    EXPECT_EQ(time.time_index(), 11);
    EXPECT_NEAR(time.new_time(), 5.5_rt, file_close(5.5_rt));

    const auto rows = read_diagnostics();
    ASSERT_EQ(rows.size(), 10);

    // Step 2: one sample, the window is filling
    EXPECT_NEAR(rows[0].time, 1.0_rt, file_close(1.0_rt));
    EXPECT_EQ(rows[0].samples, 1);
    EXPECT_EQ(rows[0].window_full, 0);
    EXPECT_EQ(rows[0].num_points, 2);

    // Step 5: the window is full on the first plateau, both points pass and
    // the hold starts. The worst speed spread is the v wobble alone, which
    // pins the sampled components to u and v
    EXPECT_NEAR(rows[3].time, 2.5_rt, file_close(2.5_rt));
    EXPECT_EQ(rows[3].samples, 4);
    EXPECT_EQ(rows[3].window_full, 1);
    EXPECT_EQ(rows[3].num_converged, 2);
    EXPECT_NEAR(
        rows[3].worst_vel_spread, scripted_speed_spread(2, 5),
        file_close(scripted_speed_spread(2, 5)));
    EXPECT_NEAR(rows[3].worst_vel_tol, 0.1_rt, file_close(0.1_rt));
    EXPECT_NEAR(rows[3].hold_elapsed, 0.0_rt, file_close(0.0_rt));

    // Step 6: the drop enters the window, the test fails and the hold resets
    EXPECT_NEAR(rows[4].time, 3.0_rt, file_close(3.0_rt));
    EXPECT_EQ(rows[4].num_converged, 0);
    EXPECT_NEAR(
        rows[4].worst_vel_spread, scripted_speed_spread(3, 6),
        file_close(scripted_speed_spread(3, 6)));
    EXPECT_NEAR(rows[4].hold_elapsed, 0.0_rt, file_close(0.0_rt));

    // Step 9: the second plateau fills the window, the hold starts again
    EXPECT_NEAR(rows[7].time, 4.5_rt, file_close(4.5_rt));
    EXPECT_EQ(rows[7].num_converged, 2);
    EXPECT_NEAR(rows[7].hold_elapsed, 0.0_rt, file_close(0.0_rt));

    // Step 11: the hold is complete. The tke spread is the scripted wobble,
    // which pins the tke sample to the tke field
    EXPECT_NEAR(rows[9].time, 5.5_rt, file_close(5.5_rt));
    EXPECT_EQ(rows[9].num_converged, 2);
    EXPECT_NEAR(rows[9].hold_elapsed, 1.0_rt, file_close(1.0_rt));
    EXPECT_NEAR(rows[9].worst_tke_spread, 0.002_rt, file_close(0.002_rt));
    EXPECT_NEAR(rows[9].worst_tke_tol, 0.01_rt, file_close(0.01_rt));
}

TEST_F(RANSConvergenceMeshTest, report_only_keeps_running)
{
    write_probe_file(100.0_rt, 300.0_rt);
    m_stop_on_convergence = false;
    setup_sim();

    auto& post_manager = sim().post_manager();
    auto& time = sim().time();
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.5_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        const auto sv = scripted_values(time.time_index());
        set_fields(sv.u, sv.v, sv.w, sv.k);
        post_manager.post_advance_work();
    }

    // The run goes to max_step, and the hold keeps counting past the point
    // where a stopping monitor would have ended it
    EXPECT_FALSE(time.stop_requested());
    EXPECT_EQ(time.time_index(), 14);

    const auto rows = read_diagnostics();
    ASSERT_EQ(rows.size(), 13);
    EXPECT_NEAR(rows[9].hold_elapsed, 1.0_rt, file_close(1.0_rt));
    EXPECT_NEAR(rows[12].hold_elapsed, 2.5_rt, file_close(2.5_rt));
}

TEST_F(RANSConvergenceMeshTest, window_is_never_full_below_min_samples)
{
    // A 1.5 s window sampled every 0.5 s holds four samples at most, so a
    // minimum of six can never be met: the window must never be reported
    // full and the run must go to max_step, however steady the history
    write_probe_file(100.0_rt, 300.0_rt);
    m_min_samples = 6;
    setup_sim();

    auto& post_manager = sim().post_manager();
    auto& time = sim().time();
    post_manager.pre_init_actions();
    post_manager.post_init_actions();

    while (time.new_timestep()) {
        time.set_current_cfl(0.45_rt / 0.5_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        const auto sv = scripted_values(time.time_index());
        set_fields(sv.u, sv.v, sv.w, sv.k);
        post_manager.post_advance_work();
    }

    EXPECT_FALSE(time.stop_requested());
    EXPECT_EQ(time.time_index(), 14);

    const auto rows = read_diagnostics();
    ASSERT_EQ(rows.size(), 13);
    for (const auto& r : rows) {
        EXPECT_EQ(r.window_full, 0);
        EXPECT_LE(r.samples, 4);
        EXPECT_NEAR(r.hold_elapsed, 0.0_rt, file_close(0.0_rt));
    }
}

TEST_F(RANSConvergenceMeshTest, rejects_a_point_below_the_terrain)
{
    // The lower point sits below a flat terrain surface at z = 150
    write_probe_file(100.0_rt, 300.0_rt);
    setup_sim();
    auto& terrain = sim().repo().declare_field("terrain_height", 1, 1);
    terrain.setVal(150.0_rt);

    auto& post_manager = sim().post_manager();
    post_manager.pre_init_actions();
    EXPECT_THROW(post_manager.post_init_actions(), std::runtime_error);
}

TEST_F(RANSConvergenceMeshTest, accepts_points_above_the_terrain)
{
    write_probe_file(200.0_rt, 300.0_rt);
    setup_sim();
    auto& terrain = sim().repo().declare_field("terrain_height", 1, 1);
    terrain.setVal(150.0_rt);

    auto& post_manager = sim().post_manager();
    post_manager.pre_init_actions();
    EXPECT_NO_THROW(post_manager.post_init_actions());
}

TEST_F(RANSConvergenceMeshTest, rejects_a_window_no_wider_than_the_interval)
{
    write_probe_file(100.0_rt, 300.0_rt);
    // Equal to the sample interval, so the window would hold one sample
    m_window = 0.5_rt;
    setup_sim();

    auto& post_manager = sim().post_manager();
    post_manager.pre_init_actions();
    EXPECT_THROW(post_manager.post_init_actions(), std::runtime_error);
}

TEST_F(RANSConvergenceMeshTest, rejects_a_zero_absolute_tolerance)
{
    write_probe_file(100.0_rt, 300.0_rt);
    m_tke_abs_tol = 0.0_rt;
    setup_sim();

    auto& post_manager = sim().post_manager();
    post_manager.pre_init_actions();
    EXPECT_THROW(post_manager.post_init_actions(), std::runtime_error);
}

TEST_F(RANSConvergenceMeshTest, rejects_a_turbulence_model_other_than_klaxell)
{
    write_probe_file(100.0_rt, 300.0_rt);
    m_turbulence_model = "Smagorinsky";
    setup_sim();

    auto& post_manager = sim().post_manager();
    post_manager.pre_init_actions();
    EXPECT_THROW(post_manager.post_init_actions(), std::runtime_error);
}

} // namespace kynema_sgf_tests
