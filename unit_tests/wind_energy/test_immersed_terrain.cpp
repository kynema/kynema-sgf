#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/iter_tools.H"
#include "ks_test_utils/test_utils.H"
#include "src/physics/ImmersedTerrain.H"
#include "AMReX_REAL.H"

#include <fstream>
#include <limits>
#include <string>

using namespace amrex::literals;

namespace {
// 100 m plateau for x in [449, 576], flat ground elsewhere
void write_terrain(const std::string& fname)
{
    std::ofstream os(fname);
    os << "6\n";
    os << "2\n";
    os << "0.0\n";
    os << "448.0\n";
    os << "449.0\n";
    os << "576.0\n";
    os << "577.0\n";
    os << "1024.0\n";
    os << "0.0\n";
    os << "1024.0\n";
    os << "0.0\n";
    os << "0.0\n";
    os << "0.0\n";
    os << "0.0\n";
    os << "100.0\n";
    os << "100.0\n";
    os << "100.0\n";
    os << "100.0\n";
    os << "0.0\n";
    os << "0.0\n";
    os << "0.0\n";
    os << "0.0\n";
}

} // namespace

namespace kynema_sgf_tests {

// Partial terrain fraction, mask classification and surface height
class ImmersedTerrainTest : public MeshTest
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
    std::string m_terrain_fname = "terrain.amrwind";
};

TEST_F(ImmersedTerrainTest, fraction_mask_and_height)
{
    using Terrain = kynema_sgf::immersedterrain::ImmersedTerrain;
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    write_terrain(m_terrain_fname);
    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();
    Terrain terrain(sim());
    const int nlevels = sim().repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = sim().repo().mesh().Geom(lev);
        terrain.initialize_fields(lev, geom);
    }
    const auto& fraction = sim().repo().get_field("terrain_fraction");
    const auto& mask = sim().repo().get_int_field("terrain_mask");
    const auto& surface = sim().repo().get_field("terrain_surface");

    // dx = dy = dz = 32. Plateau covers cells i = 14..17; height 100 fills
    // k = 0..2 and 4/32 of k = 3.

    // Fluid far from terrain
    EXPECT_NEAR(utils::field_probe(fraction, 0, 5, 5, 1), 0.0_rt, tol);
    EXPECT_EQ(utils::field_probe(mask, 0, 5, 5, 1), Terrain::mask_fluid);

    // Solid interior of the plateau
    EXPECT_NEAR(utils::field_probe(fraction, 0, 15, 10, 1), 1.0_rt, tol);
    EXPECT_EQ(utils::field_probe(mask, 0, 15, 10, 1), Terrain::mask_solid);

    // Partial cell on top of the plateau: 4 m of a 32 m cell
    EXPECT_NEAR(utils::field_probe(fraction, 0, 15, 10, 3), 0.125_rt, tol);
    EXPECT_EQ(utils::field_probe(mask, 0, 15, 10, 3), Terrain::mask_surface);

    // Cell above the partial cell is plain fluid
    EXPECT_NEAR(utils::field_probe(fraction, 0, 15, 10, 4), 0.0_rt, tol);
    EXPECT_EQ(utils::field_probe(mask, 0, 15, 10, 4), Terrain::mask_fluid);

    // Fluid cell beside the vertical wall of the plateau is a surface cell,
    // found only by the six-way neighbor search
    EXPECT_NEAR(utils::field_probe(fraction, 0, 13, 10, 1), 0.0_rt, tol);
    EXPECT_EQ(utils::field_probe(mask, 0, 13, 10, 1), Terrain::mask_surface);
    EXPECT_EQ(utils::field_probe(mask, 0, 18, 10, 1), Terrain::mask_surface);

    // Surface height on top of the plateau
    EXPECT_NEAR(
        utils::field_probe(surface, 0, 15, 10, 1, Terrain::surf_height),
        100.0_rt, tol);
}

} // namespace kynema_sgf_tests
