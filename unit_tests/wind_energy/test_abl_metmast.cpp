#include <cmath>
#include <fstream>
#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/test_utils.H"
#include "src/equation_systems/icns/icns.H"
#include "src/equation_systems/icns/source_terms/MetMastForcing.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

void write_file(const std::string& fname, const std::string& content)
{
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(fname);
        os << content;
    }
    amrex::ParallelDescriptor::Barrier();
}

// Relaxation expected at a normalized squared distance ri2 from a station
amrex::Real expected_source(
    const amrex::Real ri2,
    const amrex::Real target,
    const amrex::Real vel,
    const amrex::Real tau)
{
    return std::exp(-0.25_rt * ri2) / tau * (target - vel);
}

} // namespace

// Cells are 50 m wide and 25 m tall; centers at x, y = 25 + 50 i and
// z = 12.5 + 25 k
class MetMastTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{16, 16, 16}};
            pp.addarr("n_cell", ncell);
            pp.add("max_grid_size", 8);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> probhi{{800.0_rt, 800.0_rt, 400.0_rt}};
            pp.addarr("prob_hi", probhi);
        }
    }

    void setup_sim()
    {
        populate_parameters();
        initialize_mesh();
        auto& pde_mgr = sim().pde_manager();
        pde_mgr.register_icns();
        sim().init_physics();
        sim().repo().get_field("velocity").setVal(m_vel);
    }

    kynema_sgf::Field& src_term()
    {
        return sim().pde_manager().icns().fields().src_term;
    }

    void evaluate(const kynema_sgf::pde::icns::MetMastForcing& forcing)
    {
        src_term().setVal(0.0_rt);
        forcing(0, kynema_sgf::FieldState::New, src_term()(0));
    }

    amrex::Real probe(const int i, const int j, const int k, const int n)
    {
        return utils::field_probe(src_term(), 0, i, j, k, n);
    }

    const amrex::Vector<amrex::Real> m_vel{5.0_rt, 0.0_rt, 0.0_rt};
    const amrex::Real m_tau{30.0_rt};
    const amrex::Real m_tol{
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt};
};

TEST_F(MetMastTest, single_point)
{
    write_file("metmast_single.txt", "425 425 112.5 10 2 1 300\n");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_1dprofile_file", std::string("metmast_single.txt"));
    }
    setup_sim();
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    EXPECT_EQ(forcing.num_stations(), 1);
    evaluate(forcing);

    // At the point
    EXPECT_NEAR(probe(8, 8, 4, 0), expected_source(0, 10, 5, m_tau), m_tol);
    EXPECT_NEAR(probe(8, 8, 4, 1), expected_source(0, 2, 0, m_tau), m_tol);
    EXPECT_NEAR(probe(8, 8, 4, 2), expected_source(0, 1, 0, m_tau), m_tol);
    // Anisotropic Gaussian: 100 m away horizontally, 50 m vertically
    const amrex::Real ri2 = (0.2_rt * 0.2_rt) + (2.0_rt * 2.0_rt);
    EXPECT_NEAR(probe(10, 8, 6, 0), expected_source(ri2, 10, 5, m_tau), m_tol);
}

TEST_F(MetMastTest, multiple_points)
{
    write_file(
        "metmast_multiple.txt",
        "125 125 112.5 8 0 0 300\n675 675 112.5 2 0 0 300\n");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_1dprofile_file", std::string("metmast_multiple.txt"));
        pp.add("metmast_horizontal_radius", 50.0_rt);
    }
    setup_sim();
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    EXPECT_EQ(forcing.num_stations(), 2);
    evaluate(forcing);

    // Every point is forced, not only the first one
    EXPECT_NEAR(probe(2, 2, 4, 0), expected_source(0, 8, 5, m_tau), m_tol);
    EXPECT_NEAR(probe(13, 13, 4, 0), expected_source(0, 2, 5, m_tau), m_tol);
}

TEST_F(MetMastTest, overlapping_points)
{
    write_file(
        "metmast_overlap.txt",
        "425 425 112.5 8 0 0 300\n425 425 112.5 6 0 0 300\n");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_1dprofile_file", std::string("metmast_overlap.txt"));
    }
    setup_sim();
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    evaluate(forcing);

    // The weights are capped at one and the targets averaged
    EXPECT_NEAR(probe(8, 8, 4, 0), expected_source(0, 7, 5, m_tau), m_tol);
}

TEST_F(MetMastTest, lidar_profile)
{
    write_file(
        "lidar_profile.txt",
        "425 425\n"
        "50 4 1 0 0 0 0\n"
        "100 6 1 0 0 0 0\n"
        "200 10 1 0 0 0 0\n");
    {
        amrex::ParmParse pp("ABL");
        amrex::Vector<std::string> files{"lidar_profile.txt"};
        pp.addarr("metmast_profile_files", files);
    }
    setup_sim();
    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    EXPECT_EQ(forcing.num_stations(), 1);
    evaluate(forcing);

    // Inside the profile: linear interpolation, full weight
    EXPECT_NEAR(probe(8, 8, 3, 0), expected_source(0, 5.5, 5, m_tau), m_tol);
    EXPECT_NEAR(probe(8, 8, 5, 0), expected_source(0, 7.5, 5, m_tau), m_tol);
    EXPECT_NEAR(probe(8, 8, 5, 1), expected_source(0, 1, 0, m_tau), m_tol);
    // Below the lowest gate: bottom value, tapered over 37.5 m
    EXPECT_NEAR(
        probe(8, 8, 0, 0), expected_source(1.5_rt * 1.5_rt, 4, 5, m_tau),
        m_tol);
    // Above the highest gate: top value, tapered over 62.5 m
    EXPECT_NEAR(
        probe(8, 8, 10, 0), expected_source(2.5_rt * 2.5_rt, 10, 5, m_tau),
        m_tol);
}

TEST_F(MetMastTest, sigma_band)
{
    write_file("lidar_sigma.txt", "425 425\n87.5 7 0 0 1 0 0\n");
    {
        amrex::ParmParse pp("ABL");
        amrex::Vector<std::string> files{"lidar_sigma.txt"};
        pp.addarr("metmast_profile_files", files);
    }
    setup_sim();

    // Default band of one standard deviation: forced to the band edge 6
    {
        kynema_sgf::pde::icns::MetMastForcing forcing(sim());
        evaluate(forcing);
        EXPECT_NEAR(probe(8, 8, 3, 0), expected_source(0, 6, 5, m_tau), m_tol);
    }
    // Velocity inside a wider band: no forcing
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_sigma_factor", 3.0_rt);
        kynema_sgf::pde::icns::MetMastForcing forcing(sim());
        evaluate(forcing);
        EXPECT_NEAR(probe(8, 8, 3, 0), 0.0_rt, m_tol);
    }
    // No band: plain relaxation to the mean
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_sigma_factor", 0.0_rt);
        kynema_sgf::pde::icns::MetMastForcing forcing(sim());
        evaluate(forcing);
        EXPECT_NEAR(probe(8, 8, 3, 0), expected_source(0, 7, 5, m_tau), m_tol);
    }
}

TEST_F(MetMastTest, terrain)
{
    write_file("metmast_terrain.txt", "425 425 62.5 10 0 0 300\n");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_1dprofile_file", std::string("metmast_terrain.txt"));
    }
    setup_sim();
    // Flat terrain 100 m high: the four lowest cells are inside it
    auto& terrain_height = sim().repo().declare_field("terrain_height", 1, 1);
    auto& terrain_blank = sim().repo().declare_int_field("terrain_blank", 1, 1);
    terrain_height.setVal(100.0_rt);
    terrain_blank.setVal(0);
    terrain_blank(0).setVal(
        1, amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(15, 15, 3)), 1);

    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    evaluate(forcing);

    // The point is 62.5 m above the terrain, at z = 162.5 m
    EXPECT_NEAR(probe(8, 8, 6, 0), expected_source(0, 10, 5, m_tau), m_tol);
    // First cell above the terrain is 12.5 m above ground
    EXPECT_NEAR(probe(8, 8, 4, 0), expected_source(4, 10, 5, m_tau), m_tol);
    // Cells inside the terrain are not forced
    for (int k = 0; k < 4; ++k) {
        EXPECT_NEAR(probe(8, 8, k, 0), 0.0_rt, m_tol);
    }
}

TEST_F(MetMastTest, rate_limit)
{
    write_file("metmast_rate.txt", "425 425 112.5 10 0 0 300\n");
    {
        amrex::ParmParse pp("ABL");
        pp.add("metmast_1dprofile_file", std::string("metmast_rate.txt"));
        pp.add("meso_timescale", 1.0e-2_rt);
    }
    setup_sim();
    auto& time = sim().time();
    time.new_timestep();
    time.set_current_cfl(2.0_rt, 0.0_rt, 0.0_rt);
    EXPECT_NEAR(time.delta_t(), 0.1_rt, m_tol);

    kynema_sgf::pde::icns::MetMastForcing forcing(sim());
    evaluate(forcing);

    // tau < dt: the rate is limited to 1/dt
    EXPECT_NEAR(
        probe(8, 8, 4, 0), (10.0_rt - 5.0_rt) / 0.1_rt, 1.0e2_rt * m_tol);
}

} // namespace kynema_sgf_tests
