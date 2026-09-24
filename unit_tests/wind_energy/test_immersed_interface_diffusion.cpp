#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/iter_tools.H"
#include "ks_test_utils/test_utils.H"
#include "src/physics/ImmersedTerrain.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <string>

using namespace amrex::literals;

namespace {
// 100 m plateau for x in [449, 576], flat ground elsewhere
void write_terrain(const std::string& fname, const double height = 100.0)
{
    std::ofstream os(fname);
    os << std::setprecision(17);
    os << "6\n2\n";
    os << "0.0\n448.0\n449.0\n576.0\n577.0\n1024.0\n";
    os << "0.0\n1024.0\n";
    os << "0.0\n0.0\n0.0\n0.0\n";
    for (int n = 0; n < 4; ++n) {
        os << height << "\n";
    }
    os << "0.0\n0.0\n0.0\n0.0\n";
}
} // namespace

namespace kynema_sgf_tests {

// Face factors applied to the diffusion coefficients at the terrain interface
class ImmersedInterfaceDiffusionTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
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

    void setup(const std::string& weight, const double height = 100.0)
    {
        {
            amrex::ParmParse pp("ImmersedTerrain");
            pp.remove("drag_weight");
            pp.add("drag_weight", weight);
        }
        write_terrain("terrain.amrwind", height);
        populate_parameters();
        initialize_mesh();
        sim().pde_manager().register_icns();
        sim().init_physics();
        m_terrain = std::make_unique<Terrain>(sim());
        const int nlevels = sim().repo().num_active_levels();
        for (int lev = 0; lev < nlevels; ++lev) {
            m_terrain->initialize_fields(lev, sim().repo().mesh().Geom(lev));
        }
    }

    using Terrain = kynema_sgf::immersedterrain::ImmersedTerrain;
    std::unique_ptr<Terrain> m_terrain;
    const amrex::Real m_tol{
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt};
};

// dx = dy = dz = 32; plateau cells i = 14..17, solid for k = 0..2,
// partial (beta = 0.125) at k = 3.

TEST_F(ImmersedInterfaceDiffusionTest, center)
{
    setup("center");
    const auto& fx = sim().repo().get_field("terrain_diffusion_xf");
    const auto& fz = sim().repo().get_field("terrain_diffusion_zf");
    // Lateral wall at the face: dx / (dx/2) = 2
    EXPECT_NEAR(utils::field_probe(fx, 0, 14, 10, 1), 2.0_rt, m_tol);
    // Bottom face of the partial cell (15,10,3): centre 112 m, terrain
    // 100 m, so d1 = 12 m and the factor is 32/12
    EXPECT_NEAR(utils::field_probe(fz, 0, 15, 10, 3), 32.0_rt / 12.0_rt, m_tol);
    // Face between the partial cell (a fluid cell with center weighting) and
    // the fluid cell above: untouched
    EXPECT_NEAR(utils::field_probe(fz, 0, 15, 10, 4), 1.0_rt, m_tol);
    // Solid/solid and fluid/fluid faces untouched
    EXPECT_NEAR(utils::field_probe(fz, 0, 15, 10, 2), 1.0_rt, m_tol);
    EXPECT_NEAR(utils::field_probe(fx, 0, 5, 5, 8), 1.0_rt, m_tol);
}

TEST_F(ImmersedInterfaceDiffusionTest, fraction)
{
    setup("fraction");
    const auto& fx = sim().repo().get_field("terrain_diffusion_xf");
    const auto& fz = sim().repo().get_field("terrain_diffusion_zf");
    // Lateral wall next to an entirely solid cell: dx / (dx/2) = 2
    EXPECT_NEAR(utils::field_probe(fx, 0, 14, 10, 1), 2.0_rt, m_tol);
    // Bottom face of the partial cell (15,10,3): its fluid part spans
    // 100-128 m, centroid 114 m, so d1 = 14 m and the factor is 32/14
    EXPECT_NEAR(utils::field_probe(fz, 0, 15, 10, 3), 32.0_rt / 14.0_rt, m_tol);
    // Face above the partial cell: centroid 114 m to centre 144 m, 30 m,
    // factor 32/30 = 2 / (2 - beta)
    EXPECT_NEAR(utils::field_probe(fz, 0, 15, 10, 4), 32.0_rt / 30.0_rt, m_tol);
    // Solid/solid and fluid/fluid faces untouched
    EXPECT_NEAR(utils::field_probe(fz, 0, 15, 10, 2), 1.0_rt, m_tol);
    EXPECT_NEAR(utils::field_probe(fx, 0, 5, 5, 8), 1.0_rt, m_tol);
}

// The wall flux uses the true wall distance also when the terrain surface is
// close to a cell center (center) or to the top of a cell (fraction)

TEST_F(ImmersedInterfaceDiffusionTest, center_wall_near_cell_center)
{
    // Terrain 1.6 m (0.05 dz) below the center of cell k = 3 (112 m)
    setup("center", 110.4);
    const auto& fz = sim().repo().get_field("terrain_diffusion_zf");
    const amrex::Real expected = 32.0_rt / 1.6_rt;
    EXPECT_NEAR(
        utils::field_probe(fz, 0, 15, 10, 3), expected, 2.0e-4_rt * expected);
}

TEST_F(ImmersedInterfaceDiffusionTest, fraction_wall_near_cell_top)
{
    // Partial cell k = 3 with beta = 0.985: fluid part 127.52-128 m, so
    // d1 = 0.24 m
    setup("fraction", 127.52);
    const auto& fz = sim().repo().get_field("terrain_diffusion_zf");
    const amrex::Real expected = 32.0_rt / 0.24_rt;
    EXPECT_NEAR(
        utils::field_probe(fz, 0, 15, 10, 3), expected, 2.0e-4_rt * expected);
}

} // namespace kynema_sgf_tests
