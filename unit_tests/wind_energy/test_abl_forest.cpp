#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/physics/ForestDrag.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/utilities/sampling/FieldNorms.H"
#include "AMReX_REAL.H"

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

} // namespace

namespace amr_wind_tests {

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
    amr_wind::forestdrag::ForestDrag forest_drag(sim());
    const int nlevels = sim().repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = sim().repo().mesh().Geom(lev);
        forest_drag.initialize_fields(lev, geom);
    }

    constexpr amrex::Real n_forests = 3.0_rt;
    const auto& f_id = sim().repo().get_field("forest_id");
    const amrex::Real max_id = amr_wind::field_ops::global_max_magnitude(f_id);
    EXPECT_EQ(max_id, n_forests);

    constexpr amrex::Real expected_max_drag = 0.050285714285714288_rt;
    const auto& f_drag = sim().repo().get_field("forest_drag");
    const amrex::Real max_drag =
        amr_wind::field_ops::global_max_magnitude(f_drag);
    EXPECT_NEAR(max_drag, expected_max_drag, amr_wind::constants::TIGHT_TOL);

    constexpr amrex::Real expected_norm_drag = 0.0030635155406915832_rt;
    const auto norm_drag =
        amr_wind::field_norms::FieldNorms::get_norm(f_drag, 0, 1, 2, false);
    EXPECT_NEAR(norm_drag, expected_norm_drag, amr_wind::constants::TIGHT_TOL);
}

} // namespace amr_wind_tests
