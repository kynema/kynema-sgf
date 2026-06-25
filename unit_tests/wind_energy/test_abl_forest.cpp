#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/iter_tools.H"
#include "ks_test_utils/test_utils.H"
#include "src/physics/ForestDrag.H"
#include "src/core/field_ops.H"
#include "src/utilities/output_quantities/FieldNorms.H"
#include "AMReX_REAL.H"
#include <cmath>
#include <fstream>

using namespace amrex::literals;

namespace {
void write_forest(const std::string& fname)
{
    std::ofstream os(fname);
    //! Write forest
    os << "1  512 256 45 200  0.2 6 0.8 \n";
    os << "1  512 512 35 200  0.2 6 0.8 \n";
    os << "1  512 612 75 200  0.2 6 0.8 \n";
    os << "2  512 762 120 200 0.2 10 0.8 \n";
}

void write_point_cloud_forest(const std::string& fname)
{
    std::ofstream os(fname);
    //! x y z lad
    os << "2.5 2.5 2.5 1.0\n";
    os << "4.5 2.5 2.5 3.0\n";
    // Far point with very large LAD used to confirm nearest-point selection.
    os << "2.5 6.5 2.5 100.0\n";
}

amrex::Real idw_lad_from_two(
    const amrex::Real d1,
    const amrex::Real lad1,
    const amrex::Real d2,
    const amrex::Real lad2,
    const amrex::Real eps)
{
    const amrex::Real w1 = 1.0_rt / std::sqrt((d1 * d1) + (eps * eps));
    const amrex::Real w2 = 1.0_rt / std::sqrt((d2 * d2) + (eps * eps));
    return ((w1 * lad1) + (w2 * lad2)) / (w1 + w2);
}

} // namespace

namespace kynema_sgf_tests {

// Testing the terrain drag reading to ensure that terrain is properly setup
class ForestTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        // Make computational domain like ABL mesh
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 32, 16}};
            pp.addarr("n_cell", ncell);
            pp.add("blocking_factor", 2);
        }

        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> probhi{{1024.0_rt, 1024.0_rt, 512.0_rt}};
            pp.addarr("prob_hi", probhi);
        }
    }
    std::string m_forest_fname = "forest.amrwind";
};

class PointCloudForestTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{8, 8, 8}};
            pp.addarr("n_cell", ncell);
            pp.add("blocking_factor", 2);
        }

        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> probhi{{8.0_rt, 8.0_rt, 8.0_rt}};
            pp.addarr("prob_hi", probhi);
        }
    }

    std::string m_point_cloud_fname{"forest_points.dat"};
};

TEST_F(ForestTest, forest)
{
    // Write target wind file
    write_forest(m_forest_fname);
    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();
    amrex::ParmParse pp("incflo");
    amrex::Vector<std::string> physics{"forestDrag"};
    pp.addarr("physics", physics);
    kynema_sgf::forestdrag::ForestDrag forest_drag(sim());
    const int nlevels = sim().repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = sim().repo().mesh().Geom(lev);
        forest_drag.initialize_fields(lev, geom);
    }

    constexpr amrex::Real n_forests = 3.0_rt;
    const auto& f_id = sim().repo().get_field("forest_id");
    const amrex::Real max_id =
        kynema_sgf::field_ops::global_max_magnitude(f_id);
    EXPECT_EQ(max_id, n_forests);

    constexpr amrex::Real expected_max_drag = 0.050285714285714288_rt;
    const auto& f_drag = sim().repo().get_field("forest_drag");
    const amrex::Real max_drag =
        kynema_sgf::field_ops::global_max_magnitude(f_drag);
    EXPECT_NEAR(max_drag, expected_max_drag, kynema_sgf::constants::TIGHT_TOL);

    constexpr amrex::Real expected_norm_drag = 0.0030635155406915832_rt;
    const auto norm_drag =
        kynema_sgf::field_norms::FieldNorms::get_norm(f_drag, 0, 1, 2, false);
    EXPECT_NEAR(
        norm_drag, expected_norm_drag, kynema_sgf::constants::TIGHT_TOL);
}

TEST_F(PointCloudForestTest, point_cloud_selection_and_interpolation)
{
    write_point_cloud_forest(m_point_cloud_fname);
    populate_parameters();
    initialize_mesh();

    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    amrex::ParmParse pp("ForestDrag");
    amrex::Vector<std::string> cloud_files{m_point_cloud_fname};
    amrex::Vector<amrex::Real> cds{2.0_rt};
    const amrex::Real tol = kynema_sgf::constants::TIGHT_TOL;
    pp.addarr("point_cloud_files", cloud_files);
    pp.addarr("coefficients_of_drag", cds);
    pp.add("point_neighbors", 2);
    pp.add("point_interp_eps", tol);

    kynema_sgf::forestdrag::ForestDrag forest_drag(sim());
    forest_drag.initialize_fields(0, sim().repo().mesh().Geom(0));

    const auto& f_drag = sim().repo().get_field("forest_drag");
    const auto& f_id = sim().repo().get_field("forest_id");

    // Exact-point selections: drag = cd * lad
    // Point cloud forms convex hull triangle in xy-plane:
    // Vertices: (2.5, 2.5), (4.5, 2.5), (2.5, 6.5)
    const auto drag_exact_p1 =
        utils::field_probe(f_drag, 0, 2, 2, 2); // x,y,z = 2.5,2.5,2.5
    const auto drag_exact_p2 =
        utils::field_probe(f_drag, 0, 4, 2, 2); // x,y,z = 4.5,2.5,2.5
    EXPECT_NEAR(drag_exact_p1, 2.0_rt, tol);
    EXPECT_NEAR(drag_exact_p2, 6.0_rt, tol);
    EXPECT_NEAR(utils::field_probe(f_id, 0, 2, 2, 2), 0.0_rt, tol);

    // Midpoint interpolation from the two nearest points (far LAD=100 point
    // must be ignored because neighbors=2). Cell is inside hull.
    const auto drag_mid =
        utils::field_probe(f_drag, 0, 3, 2, 2); // x,y,z = 3.5,2.5,2.5
    EXPECT_NEAR(drag_mid, 4.0_rt, tol);

    // Unequal-distance interpolation (still only nearest two points).
    // Cell at (2.5, 3.5) is inside hull.
    const auto drag_off =
        utils::field_probe(f_drag, 0, 2, 3, 2); // x,y,z = 2.5,3.5,2.5
    const auto expected_lad =
        idw_lad_from_two(1.0_rt, 1.0_rt, std::sqrt(5.0_rt), 3.0_rt, tol);
    EXPECT_NEAR(drag_off, 2.0_rt * expected_lad, tol);

    // Cell at (3.5, 3.5) should be inside hull and interpolate.
    const auto drag_interior =
        utils::field_probe(f_drag, 0, 3, 3, 2); // x,y,z = 3.5,3.5,2.5
    EXPECT_GT(
        drag_interior, 0.0_rt); // Should have non-zero drag if inside hull
    EXPECT_NEAR(utils::field_probe(f_id, 0, 3, 3, 2), 0.0_rt, tol);

    // Outside all cloud extents should remain untouched.
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 0, 0, 0), 0.0_rt, tol);
    EXPECT_NEAR(utils::field_probe(f_id, 0, 0, 0, 0), -1.0_rt, tol);

    // Cell at (5, 2) center (5.5, 2.5) is outside hull (x > 4.5).
    // Even though it's near the second cloud point, it should be excluded.
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 5, 2, 2), 0.0_rt, tol);
    EXPECT_NEAR(utils::field_probe(f_id, 0, 5, 2, 2), -1.0_rt, tol);

    // Cell at (0, 2) center (0.5, 2.5) is outside hull (x < 2.5).
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 0, 2, 2), 0.0_rt, tol);
    EXPECT_NEAR(utils::field_probe(f_id, 0, 0, 2, 2), -1.0_rt, tol);

    // Cell at (2, 7) center (2.5, 7.5) is outside hull (y > 6.5).
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 2, 7, 2), 0.0_rt, tol);
    EXPECT_NEAR(utils::field_probe(f_id, 0, 2, 7, 2), -1.0_rt, tol);

    // Vertical canopy cutoff: all cloud points are at z=2.5, so
    // max_z_neighbors=2.5. The guard is (z - 0.5*dz) <= max_z_neighbors.
    // dx[2]=1.0, so k=2 -> cell bottom 2.0 <= 2.5 (inside), k=3 -> 3.0 > 2.5
    // (above canopy). Positive control: cell (3,2,2) center (3.5,2.5,2.5) is
    // inside hull and within canopy height — must have non-zero drag.
    EXPECT_GT(utils::field_probe(f_drag, 0, 3, 2, 2), 0.0_rt);

    // Above-canopy cell (3,2,3) center (3.5,2.5,3.5): x-y is inside hull but
    // the vertical check zeroes the contribution because
    // z - 0.5*dz = 3.0 > max_z_neighbors = 2.5.
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 3, 2, 3), 0.0_rt, tol);
    EXPECT_NEAR(utils::field_probe(f_id, 0, 3, 2, 3), -1.0_rt, tol);
}

} // namespace kynema_sgf_tests
