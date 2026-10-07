#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/iter_tools.H"
#include "ks_test_utils/test_utils.H"
#include "src/physics/ForestDrag.H"
#include "src/physics/TerrainDrag.H"
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

// One type-2 forest whose top (48 m) is a cell center of the 32 m grid
void write_type2_forest_top_on_center(const std::string& fname)
{
    std::ofstream os(fname);
    os << "2  512 512 48 200 0.2 6 0.8 \n";
}

// Terrain on a 2 x 2 grid, linear in x: height = z_west + slope * (x - x_west)
void write_linear_terrain(
    const std::string& fname,
    const amrex::Real x_west,
    const amrex::Real x_east,
    const amrex::Real y_south,
    const amrex::Real y_north,
    const amrex::Real z_west,
    const amrex::Real z_east)
{
    std::ofstream os(fname);
    os << "2\n2\n";
    os << x_west << "\n" << x_east << "\n";
    os << y_south << "\n" << y_north << "\n";
    // z index is (i * ny) + j
    os << z_west << "\n" << z_west << "\n";
    os << z_east << "\n" << z_east << "\n";
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
    std::string m_terrain_fname = "forest_terrain.amrwind";
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
    std::string m_terrain_fname{"forest_cloud_terrain.amrwind"};
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

// Flat terrain 128 m (4 cells) high: every forest moves up by 4 cells and
// keeps the legacy maximum and norm.
TEST_F(ForestTest, forest_on_flat_terrain)
{
    write_forest(m_forest_fname);
    write_linear_terrain(
        m_terrain_fname, 0.0_rt, 1024.0_rt, 0.0_rt, 1024.0_rt, 128.0_rt,
        128.0_rt);
    populate_parameters();
    {
        amrex::ParmParse pp("TerrainDrag");
        pp.add("terrain_file", m_terrain_fname);
    }
    initialize_mesh();
    sim().pde_manager().register_icns();
    auto& terrain_drag = sim().physics_manager().create("TerrainDrag", sim());
    auto& forest_drag = sim().physics_manager().create("ForestDrag", sim());
    const auto& geom = sim().repo().mesh().Geom(0);
    terrain_drag.initialize_fields(0, geom);
    forest_drag.initialize_fields(0, geom);

    const amrex::Real tol = kynema_sgf::constants::TIGHT_TOL;
    const auto& f_id = sim().repo().get_field("forest_id");
    const auto& f_drag = sim().repo().get_field("forest_drag");
    EXPECT_EQ(kynema_sgf::field_ops::global_max_magnitude(f_id), 3.0_rt);
    EXPECT_NEAR(
        kynema_sgf::field_ops::global_max_magnitude(f_drag),
        0.050285714285714288_rt, tol);
    EXPECT_NEAR(
        kynema_sgf::field_norms::FieldNorms::get_norm(f_drag, 0, 1, 2, false),
        0.0030635155406915832_rt, tol);

    // Forest 0 (type 1, height 45 m, cd 0.2, lai 6) around x = 512, y = 256.
    // Cell centers z = 16 + 32 k; ground at 128 m, so the canopy fills
    // k = 4 (16 m) and k = 5 (48 m > 45 m is above it).
    const amrex::Real lad = 0.2_rt * 6.0_rt / 45.0_rt;
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 16, 8, 0), 0.0_rt, tol);
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 16, 8, 3), 0.0_rt, tol);
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 16, 8, 4), lad, tol);
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 16, 8, 5), 0.0_rt, tol);
    EXPECT_NEAR(utils::field_probe(f_id, 0, 16, 8, 3), -1.0_rt, tol);
    EXPECT_NEAR(utils::field_probe(f_id, 0, 16, 8, 4), 0.0_rt, tol);
}

// Terrain rising 0.25 m per m in x: the canopy base follows the local ground
// in each column.
TEST_F(ForestTest, forest_on_sloped_terrain)
{
    write_forest(m_forest_fname);
    write_linear_terrain(
        m_terrain_fname, 0.0_rt, 1024.0_rt, 0.0_rt, 1024.0_rt, 0.0_rt,
        256.0_rt);
    populate_parameters();
    {
        amrex::ParmParse pp("TerrainDrag");
        pp.add("terrain_file", m_terrain_fname);
    }
    initialize_mesh();
    sim().pde_manager().register_icns();
    auto& terrain_drag = sim().physics_manager().create("TerrainDrag", sim());
    auto& forest_drag = sim().physics_manager().create("ForestDrag", sim());
    const auto& geom = sim().repo().mesh().Geom(0);
    terrain_drag.initialize_fields(0, geom);
    forest_drag.initialize_fields(0, geom);

    const amrex::Real tol = kynema_sgf::constants::TIGHT_TOL;
    const auto& f_drag = sim().repo().get_field("forest_drag");
    const amrex::Real lad = 0.2_rt * 6.0_rt / 45.0_rt;

    // i = 13: x = 432, ground 108 m, canopy 108-153 m -> k = 3, 4.
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 13, 8, 2), 0.0_rt, tol);
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 13, 8, 3), lad, tol);
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 13, 8, 4), lad, tol);
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 13, 8, 5), 0.0_rt, tol);

    // i = 16: x = 528, ground 132 m, canopy 132-177 m -> k = 4, 5.
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 16, 8, 3), 0.0_rt, tol);
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 16, 8, 4), lad, tol);
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 16, 8, 5), lad, tol);
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 16, 8, 6), 0.0_rt, tol);

    // The tallest forest (120 m) on the highest ground still fits in the
    // bounding box, which is extended by the terrain maximum.
    EXPECT_EQ(
        kynema_sgf::field_ops::global_max_magnitude(
            sim().repo().get_field("forest_id")),
        3.0_rt);
}

// ForestDrag listed before TerrainDrag would read an unset terrain_height
TEST_F(ForestTest, forest_before_terrain_aborts)
{
    write_forest(m_forest_fname);
    write_linear_terrain(
        m_terrain_fname, 0.0_rt, 1024.0_rt, 0.0_rt, 1024.0_rt, 128.0_rt,
        128.0_rt);
    populate_parameters();
    {
        amrex::ParmParse pp("TerrainDrag");
        pp.add("terrain_file", m_terrain_fname);
    }
    initialize_mesh();
    sim().pde_manager().register_icns();
    auto& forest_drag = sim().physics_manager().create("ForestDrag", sim());
    auto& terrain_drag = sim().physics_manager().create("TerrainDrag", sim());
    const auto& geom = sim().repo().mesh().Geom(0);
    terrain_drag.initialize_fields(0, geom);
    EXPECT_THROW(forest_drag.initialize_fields(0, geom), amrex::RuntimeError);
}

// terrain_aware = false keeps the legacy placement with TerrainDrag active
TEST_F(ForestTest, forest_terrain_aware_off)
{
    write_forest(m_forest_fname);
    write_linear_terrain(
        m_terrain_fname, 0.0_rt, 1024.0_rt, 0.0_rt, 1024.0_rt, 128.0_rt,
        128.0_rt);
    populate_parameters();
    {
        amrex::ParmParse pp("TerrainDrag");
        pp.add("terrain_file", m_terrain_fname);
    }
    {
        amrex::ParmParse pp("ForestDrag");
        pp.add("terrain_aware", false);
    }
    initialize_mesh();
    sim().pde_manager().register_icns();
    auto& terrain_drag = sim().physics_manager().create("TerrainDrag", sim());
    auto& forest_drag = sim().physics_manager().create("ForestDrag", sim());
    const auto& geom = sim().repo().mesh().Geom(0);
    terrain_drag.initialize_fields(0, geom);
    forest_drag.initialize_fields(0, geom);

    const amrex::Real tol = kynema_sgf::constants::TIGHT_TOL;
    const auto& f_drag = sim().repo().get_field("forest_drag");
    EXPECT_NEAR(
        utils::field_probe(f_drag, 0, 16, 8, 0), 0.2_rt * 6.0_rt / 45.0_rt,
        tol);
    EXPECT_NEAR(
        kynema_sgf::field_norms::FieldNorms::get_norm(f_drag, 0, 1, 2, false),
        0.0030635155406915832_rt, tol);
}

// Terrain built from single-phase OceanWaves is filled after the levels are
// set up and moves every step, so the forests keep their legacy placement.
TEST_F(ForestTest, forest_ignores_wave_terrain)
{
    write_forest(m_forest_fname);
    populate_parameters();
    {
        amrex::ParmParse pp("OceanWaves");
        pp.add("label", (std::string) "lin_ow");
        amrex::ParmParse ppow("OceanWaves.lin_ow");
        ppow.add("type", (std::string) "LinearWaves");
        ppow.add("wave_height", 2.0_rt);
        ppow.add("wave_length", 40.0_rt);
        ppow.add("water_depth", 133.0_rt);
    }
    initialize_mesh();
    sim().pde_manager().register_icns();
    sim().physics_manager().create("OceanWaves", sim());
    auto& terrain_drag = sim().physics_manager().create("TerrainDrag", sim());
    auto& forest_drag = sim().physics_manager().create("ForestDrag", sim());
    EXPECT_TRUE(
        sim()
            .physics_manager()
            .get<kynema_sgf::terraindrag::TerrainDrag>()
            .terrain_is_waves());

    // Stand-in for a wave surface, which TerrainDrag::initialize_fields does
    // not fill
    sim().repo().get_field("terrain_height").setVal(128.0_rt);
    const auto& geom = sim().repo().mesh().Geom(0);
    terrain_drag.initialize_fields(0, geom);
    forest_drag.initialize_fields(0, geom);

    const amrex::Real tol = kynema_sgf::constants::TIGHT_TOL;
    const auto& f_drag = sim().repo().get_field("forest_drag");
    EXPECT_NEAR(
        utils::field_probe(f_drag, 0, 16, 8, 0), 0.2_rt * 6.0_rt / 45.0_rt,
        tol);
    EXPECT_NEAR(
        kynema_sgf::field_norms::FieldNorms::get_norm(f_drag, 0, 1, 2, false),
        0.0030635155406915832_rt, tol);
}

// Point-cloud heights are above the local ground: 2 m of flat terrain moves
// the samples at z = 2.5 to the cells at z = 4.5.
TEST_F(PointCloudForestTest, point_cloud_on_flat_terrain)
{
    write_point_cloud_forest(m_point_cloud_fname);
    write_linear_terrain(
        m_terrain_fname, 0.0_rt, 8.0_rt, 0.0_rt, 8.0_rt, 2.0_rt, 2.0_rt);
    populate_parameters();
    const amrex::Real tol = kynema_sgf::constants::TIGHT_TOL;
    {
        amrex::ParmParse pp("TerrainDrag");
        pp.add("terrain_file", m_terrain_fname);
    }
    {
        amrex::ParmParse pp("ForestDrag");
        amrex::Vector<std::string> cloud_files{m_point_cloud_fname};
        amrex::Vector<amrex::Real> cds{2.0_rt};
        pp.addarr("point_cloud_files", cloud_files);
        pp.addarr("coefficients_of_drag", cds);
        pp.add("point_neighbors", 2);
        pp.add("point_interp_eps", tol);
    }
    initialize_mesh();
    sim().pde_manager().register_icns();
    auto& terrain_drag = sim().physics_manager().create("TerrainDrag", sim());
    auto& forest_drag = sim().physics_manager().create("ForestDrag", sim());
    const auto& geom = sim().repo().mesh().Geom(0);
    terrain_drag.initialize_fields(0, geom);
    forest_drag.initialize_fields(0, geom);

    const auto& f_drag = sim().repo().get_field("forest_drag");
    const auto& f_id = sim().repo().get_field("forest_id");

    // Samples, shifted 2 cells up
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 2, 2, 4), 2.0_rt, tol);
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 4, 2, 4), 6.0_rt, tol);
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 3, 2, 4), 4.0_rt, tol);
    EXPECT_NEAR(utils::field_probe(f_id, 0, 2, 2, 4), 0.0_rt, tol);

    // Above the canopy top (2.5 m above ground)
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 3, 2, 5), 0.0_rt, tol);

    // Inside the terrain
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 2, 2, 1), 0.0_rt, tol);
    EXPECT_NEAR(utils::field_probe(f_id, 0, 2, 2, 1), -1.0_rt, tol);
}

// The type-2 profile divided by zero at the tree top (NaN drag); it is zero
// there, in the legacy placement and on the terrain.
TEST_F(ForestTest, type2_forest_top_on_cell_center)
{
    write_type2_forest_top_on_center(m_forest_fname);
    populate_parameters();
    initialize_mesh();
    sim().pde_manager().register_icns();
    kynema_sgf::forestdrag::ForestDrag forest_drag(sim());
    forest_drag.initialize_fields(0, sim().repo().mesh().Geom(0));

    const amrex::Real tol = kynema_sgf::constants::TIGHT_TOL;
    const auto& f_drag = sim().repo().get_field("forest_drag");
    EXPECT_TRUE(
        std::isfinite(kynema_sgf::field_ops::global_max_magnitude(f_drag)));
    // k = 0 (z = 16 m) is in the canopy, k = 1 (z = 48 m) is the tree top
    EXPECT_GT(utils::field_probe(f_drag, 0, 16, 16, 0), 0.0_rt);
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 16, 16, 1), 0.0_rt, tol);
    EXPECT_NEAR(
        utils::field_probe(sim().repo().get_field("forest_id"), 0, 16, 16, 1),
        0.0_rt, tol);
}

TEST_F(ForestTest, type2_forest_top_on_cell_center_terrain)
{
    write_type2_forest_top_on_center(m_forest_fname);
    write_linear_terrain(
        m_terrain_fname, 0.0_rt, 1024.0_rt, 0.0_rt, 1024.0_rt, 128.0_rt,
        128.0_rt);
    populate_parameters();
    {
        amrex::ParmParse pp("TerrainDrag");
        pp.add("terrain_file", m_terrain_fname);
    }
    initialize_mesh();
    sim().pde_manager().register_icns();
    auto& terrain_drag = sim().physics_manager().create("TerrainDrag", sim());
    auto& forest_drag = sim().physics_manager().create("ForestDrag", sim());
    const auto& geom = sim().repo().mesh().Geom(0);
    terrain_drag.initialize_fields(0, geom);
    forest_drag.initialize_fields(0, geom);

    const amrex::Real tol = kynema_sgf::constants::TIGHT_TOL;
    const auto& f_drag = sim().repo().get_field("forest_drag");
    EXPECT_TRUE(
        std::isfinite(kynema_sgf::field_ops::global_max_magnitude(f_drag)));
    // Ground at 128 m: k = 4 is 16 m and k = 5 is 48 m above it
    EXPECT_GT(utils::field_probe(f_drag, 0, 16, 16, 4), 0.0_rt);
    EXPECT_NEAR(utils::field_probe(f_drag, 0, 16, 16, 5), 0.0_rt, tol);
}

} // namespace kynema_sgf_tests
