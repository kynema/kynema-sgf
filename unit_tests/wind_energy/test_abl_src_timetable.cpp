#include <numbers>
#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/equation_systems/icns/icns.H"
#include "amr-wind/equation_systems/icns/icns_ops.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/equation_systems/icns/source_terms/BodyForce.H"
#include "amr-wind/equation_systems/icns/source_terms/ABLForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/GeostrophicForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/CoriolisForcing.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace {
void write_target_velocity_file(const std::string& fname)
{
    std::ofstream os(fname);
    // Write header **??**
    os << "Time\tWindSpeed\tHorzAngle\n";
    // Write time table
    os << "0.0\t8.0\t0.0\n";
    os << "0.1\t8.0\t360.0\n";
    os << "0.3\t8.0\t5.0\n";
    os << "0.5\t8.0\t12.0\n";
}
void write_body_force_file(const std::string& fname)
{
    std::ofstream os(fname);
    // Write header **??**
    os << "time\tfx\tfy\tfz\n";
    // Write time table
    os << "0.000000000000000\t0.000000000000000\t0.000000000000000\t0."
          "000000000000000\n";
    os << "0.050000000000000\t0.000000000000000\t0.000000000000000\t0."
          "000000000000000\n";
    os << "0.150000000000000\t-0.076142273451376\t3.489550989226880\t0."
          "000000000000000\n";
    os << "0.250000000000000\t-0.304424152660356\t6.972459419812656\t0."
          "000000000000000\n";
    os << "0.350000000000000\t-0.878730931046654\t11.824752890368849\t0."
          "000000000000000\n";
}
} // namespace

namespace amr_wind_tests {

using ICNSFields =
    amr_wind::pde::FieldRegOp<amr_wind::pde::ICNS, amr_wind::fvm::Godunov>;

// Tests in this file involve ABLForcing, BodyForce, and GeostrophicWind
class ABLSrcTimeTableTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        // Make computational domain like ABL mesh
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{8, 8, 64}};
            pp.addarr("n_cell", ncell);
        }

        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> probhi{{120.0_rt, 120.0_rt, 1000.0_rt}};
            pp.addarr("prob_hi", probhi);
        }

        // Parameters for initializing ABL case
        {
            amrex::ParmParse pp("ABL");
            amrex::Vector<amrex::Real> theights{
                {0.0_rt, 650.0_rt, 750.0_rt, 1000.0_rt}};
            amrex::Vector<amrex::Real> tvalues{
                {300.0_rt, 300.0_rt, 308.0_rt, 308.75_rt}};
            pp.addarr("temperature_heights", theights);
            pp.addarr("temperature_values", tvalues);
            pp.add("perturb_ref_height", 50.0_rt);
            pp.add("kappa", 0.41_rt);
            pp.add("surface_roughness_z0", 0.1_rt);
        }

        {
            // Physics
            amrex::ParmParse pp("incflo");
            amrex::Vector<std::string> phystr{"ABL"};
            pp.addarr("physics", phystr);
            pp.add("density", 1.0_rt);
        }

        // Timestep size
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", 0.1_rt);
        }

        // Transport
        {
            amrex::ParmParse pp("transport");
            pp.add("reference_temperature", 300.0_rt);
        }
    }
    std::string m_tvel_fname = "target_velocities.txt";
    std::string m_forces_fname = "abl_forces.txt";
    amrex::Real m_dt{0.1_rt};
    int m_nsteps = 5;
};

TEST_F(ABLSrcTimeTableTest, abl)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    // Write target wind file
    write_target_velocity_file(m_tvel_fname);

    // Set up simulation parameters and mesh
    populate_parameters();
    // ABL Forcing
    {
        amrex::ParmParse pp("ABLForcing");
        pp.add("abl_forcing_height", 90.0_rt);
        pp.add("velocity_timetable", m_tvel_fname);
    }
    initialize_mesh();

    // Set up PDEs and physics objects
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    // Get icns source term, which will be tested
    auto& src_term = pde_mgr.icns().fields().src_term;
    src_term.setVal(0.0_rt);
    // Source term object for ABLForcing
    amr_wind::pde::icns::ABLForcing abl_forcing(sim());
    run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term(lev).array(mfi);

        abl_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });

    // Initial velocity should be the same as target velocity, so force is 0
    amrex::Vector<amrex::Real> target_force = {0.0_rt, 0.0_rt, 0.0_rt};
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const auto min_val = utils::field_min(src_term, i);
        const auto max_val = utils::field_max(src_term, i);
        EXPECT_NEAR(min_val, target_force[i], tol);
        EXPECT_NEAR(min_val, max_val, tol);
    }

    // Change mean velocity and recalculate
    const amrex::Vector<amrex::Real> init_vel{8.0_rt, 0.0_rt, 0.0_rt};
    const amrex::Vector<amrex::Real> new_vel{5.0_rt, 3.0_rt, 0.0_rt};
    // This function actually calculates the forcing terms, applied later
    // This function also is the one that writes forces to file
    abl_forcing.set_mean_velocities(new_vel[0], new_vel[1]);
    // Doing this through ABL physics is complicated and would require
    // calc_averages in ABLStats, then pre_advance_work in ABL physics
    src_term.setVal(0.0_rt);
    run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term(lev).array(mfi);

        abl_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });
    // Velocity at hub height is (5, 3) and target is (8, 0)
    target_force[0] = (init_vel[0] - new_vel[0]) / m_dt;
    target_force[1] = (init_vel[1] - new_vel[1]) / m_dt;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const auto min_val = utils::field_min(src_term, i);
        const auto max_val = utils::field_max(src_term, i);
        EXPECT_NEAR(min_val, target_force[i], tol);
        EXPECT_NEAR(min_val, max_val, tol);
    }

    // Advance time (twice to get to interpolation interval)
    sim().time().new_timestep();
    sim().time().advance_time();
    sim().time().new_timestep();
    sim().time().advance_time();

    // Go back to original mean values on mesh
    abl_forcing.set_mean_velocities(init_vel[0], init_vel[1]);
    // Recalculate forcing and check
    src_term.setVal(0.0_rt);
    run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term(lev).array(mfi);

        abl_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });
    // Velocity at hub height is 8 at 0deg and target is 8 at 2.5deg
    target_force[0] =
        (8.0_rt *
             std::cos(std::numbers::pi_v<amrex::Real> / 180.0_rt * 2.5_rt) -
         init_vel[0]) /
        m_dt;
    target_force[1] =
        (8.0_rt *
             std::sin(std::numbers::pi_v<amrex::Real> / 180.0_rt * 2.5_rt) -
         init_vel[1]) /
        m_dt;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const auto min_val = utils::field_min(src_term, i);
        const auto max_val = utils::field_max(src_term, i);
        EXPECT_NEAR(min_val, target_force[i], tol);
        EXPECT_NEAR(min_val, max_val, tol);
    }

    // Delete target wind file
    const char* fname = m_tvel_fname.c_str();
    {
        std::ifstream f(fname);
        if (f.good()) {
            remove(fname);
        }
        // Check that file is removed
        std::ifstream ff(fname);
        EXPECT_FALSE(ff.good());
    }
}

TEST_F(ABLSrcTimeTableTest, bodyforce)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    // Write bodyforce file
    write_body_force_file(m_forces_fname);

    // Set up simulation parameters and mesh
    populate_parameters();
    // Body Forcing
    {
        amrex::ParmParse pp("BodyForce");
        pp.add("type", (std::string) "uniform_timetable");
        pp.add("uniform_timetable_file", m_forces_fname);
    }
    // Velocity
    {
        amrex::ParmParse pp("incflo");
        pp.addarr(
            "velocity", amrex::Vector<amrex::Real>{8.0_rt, 0.0_rt, 0.0_rt});
    }
    initialize_mesh();

    // Set up PDEs and physics objects
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    // Get icns source term, which will be tested
    auto& src_term = pde_mgr.icns().fields().src_term;
    src_term.setVal(0.0_rt);
    // Source term object for BodyForce
    amr_wind::pde::icns::BodyForce body_forcing(sim());
    run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term(lev).array(mfi);

        body_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });

    // Force is 0 initially
    amrex::Vector<amrex::Real> target_force = {0.0_rt, 0.0_rt, 0.0_rt};
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const auto min_val = utils::field_min(src_term, i);
        const auto max_val = utils::field_max(src_term, i);
        EXPECT_NEAR(min_val, target_force[i], tol);
        EXPECT_NEAR(min_val, max_val, tol);
    }

    amrex::Vector<amrex::Real> angles{0.0_rt, 2.5_rt, 5.0_rt, 8.5_rt};

    // Loop through timesteps to check force outputs
    for (int n = 0; n < m_nsteps - 1; ++n) {
        // Advance time
        sim().time().new_timestep();
        sim().time().advance_time();
        // Recalculate forcing and check
        src_term.setVal(0.0_rt);
        run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
            const auto& bx = mfi.tilebox();
            const auto& src_arr = src_term(lev).array(mfi);

            body_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
        });
        // Forces correspond to ABL Forcing from other test
        const amrex::Vector<amrex::Real> init_vel{8.0_rt, 0.0_rt, 0.0_rt};
        target_force[0] = (8.0_rt * std::cos(
                                        std::numbers::pi_v<amrex::Real> /
                                        180.0_rt * angles[n]) -
                           init_vel[0]) /
                          m_dt;
        target_force[1] = (8.0_rt * std::sin(
                                        std::numbers::pi_v<amrex::Real> /
                                        180.0_rt * angles[n]) -
                           init_vel[1]) /
                          m_dt;
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            const auto min_val = utils::field_min(src_term, i);
            const auto max_val = utils::field_max(src_term, i);
            EXPECT_NEAR(min_val, target_force[i], tol);
            EXPECT_NEAR(min_val, max_val, tol);
        }
    }

    // Delete body force file
    const char* fname = m_forces_fname.c_str();
    {
        std::ifstream f(fname);
        if (f.good()) {
            remove(fname);
        }
        // Check that file is removed
        std::ifstream ff(fname);
        EXPECT_FALSE(ff.good());
    }
}

TEST_F(ABLSrcTimeTableTest, geostrophic)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    // Default Coriolis parameters
    const amrex::Real cf =
        2.0_rt * (2.0_rt * std::numbers::pi_v<amrex::Real> / 86164.091_rt) *
        std::sin(std::numbers::pi_v<amrex::Real> / 180.0_rt * 90.0_rt);

    // Write target wind file
    write_target_velocity_file(m_tvel_fname);

    // Set up simulation parameters and mesh
    populate_parameters();
    // Geostrophic Forcing
    {
        amrex::ParmParse pp("GeostrophicForcing");
        pp.add("geostrophic_wind_timetable", m_tvel_fname);
    }
    // Coriolis Forcing (for parameters)
    {
        amrex::ParmParse pp("CoriolisForcing");
        pp.add("latitude", 90.0_rt);
    }
    // Velocity
    {
        amrex::ParmParse pp("incflo");
        pp.addarr(
            "velocity", amrex::Vector<amrex::Real>{8.0_rt, 0.0_rt, 0.0_rt});
    }
    initialize_mesh();

    // Set up PDEs and physics objects
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    // Get icns source term, which will be tested
    auto& src_term = pde_mgr.icns().fields().src_term;
    src_term.setVal(0.0_rt);
    // Source term object for geostrophic forcing
    amr_wind::pde::icns::GeostrophicForcing gstr_forcing(sim());
    run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term(lev).array(mfi);

        gstr_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });

    // Initial velocity should be the same as target velocity
    const amrex::Vector<amrex::Real> init_vel{8.0_rt, 0.0_rt, 0.0_rt};
    amrex::Vector<amrex::Real> target_force = {
        -cf * init_vel[1], cf * init_vel[0], 0.0_rt};
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const auto min_val = utils::field_min(src_term, i);
        const auto max_val = utils::field_max(src_term, i);
        EXPECT_NEAR(min_val, target_force[i], tol);
        EXPECT_NEAR(min_val, max_val, tol);
    }

    // Advance time (twice to get to interpolation interval)
    sim().time().new_timestep();
    sim().time().advance_time();
    sim().time().new_timestep();
    sim().time().advance_time();

    // Recalculate forcing and check
    src_term.setVal(0.0_rt);
    run_algorithm(src_term, [&](const int lev, const amrex::MFIter& mfi) {
        const auto& bx = mfi.tilebox();
        const auto& src_arr = src_term(lev).array(mfi);

        gstr_forcing(lev, mfi, bx, amr_wind::FieldState::New, src_arr);
    });
    // New target velocity is 8 at 1.25deg (nph velocity)
    const amrex::Vector<amrex::Real> targ_vel{
        8.0_rt * std::cos(std::numbers::pi_v<amrex::Real> / 180.0_rt * 1.25_rt),
        8.0_rt * std::sin(std::numbers::pi_v<amrex::Real> / 180.0_rt * 1.25_rt),
        0.0_rt};
    target_force[0] = -cf * targ_vel[1];
    target_force[1] = cf * targ_vel[0];
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        const auto min_val = utils::field_min(src_term, i);
        const auto max_val = utils::field_max(src_term, i);
        EXPECT_NEAR(min_val, target_force[i], tol);
        EXPECT_NEAR(min_val, max_val, tol);
    }

    // Delete target wind file
    const char* fname = m_tvel_fname.c_str();
    {
        std::ifstream f(fname);
        if (f.good()) {
            remove(fname);
        }
        // Check that file is removed
        std::ifstream ff(fname);
        EXPECT_FALSE(ff.good());
    }
}

} // namespace amr_wind_tests
