#include <sstream>

#include "aw_test_utils/AmrexTest.H"
#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/OutputCapture.H"

#include "AMReX_Box.H"
#include "AMReX_BoxArray.H"
#include "AMReX_BoxList.H"
#include "AMReX_Geometry.H"
#include "AMReX_RealBox.H"
#include "AMReX_Vector.H"

#include "amr-wind/utilities/tagging/CartBoxRefinement.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

//! Custom test fixture for Cartesian Box refinement
class NestRefineTest : public MeshTest
{
protected:
    void setup_refinement_inputs()
    {
        populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{16, 128, 16}};

            pp.add("max_level", 1);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{-20.0_rt, -100.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{{20.0_rt, 100.0_rt, 30.0_rt}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }

    std::stringstream m_cout_buf;
    std::streambuf* m_orig_buf{nullptr};
};

TEST_F(NestRefineTest, box_refine)
{
    setup_refinement_inputs();
    // Create the "input file"
    std::stringstream ss;
    ss << "1 // Number of levels" << '\n';
    ss << "4 // Number of boxes at this level" << '\n';
    ss << "-10.0 -75.0 0.0 15.0 -65.0 20.0" << '\n';
    ss << "-10.0 -55.0 0.0 15.0 -45.0 20.0" << '\n';
    ss << "-10.0  25.0 0.0 15.0  35.0 20.0" << '\n';
    ss << "-10.0  65.0 0.0 15.0  75.0 20.0" << '\n';

    create_mesh_instance<RefineMesh>();
    auto& ref_vec = mesh<RefineMesh>()->refine_criteria_vec();
    ref_vec.emplace_back(std::make_unique<amr_wind::CartBoxRefinement>(sim()));
    auto* box_refine =
        dynamic_cast<amr_wind::CartBoxRefinement*>(ref_vec[0].get());
    box_refine->read_inputs(mesh(), ss);
    // Store the target boxarray for future tests
    auto targets = box_refine->boxarray_vec();
    initialize_mesh();

    auto ba1 = mesh().boxArray(1);
    // NOTE: Box definitions were based on Level 0 to tag cells
    auto ba2 = targets[0].refine(2);
    EXPECT_TRUE(ba1.contains(ba2));
}

/*  Check that the implementation emits a warning when the levels requested in
 *  the input file does not match the levels in static refinement file
 */
TEST_F(NestRefineTest, level_warning)
{
    setup_refinement_inputs();
    {
        amrex::ParmParse pp("amr");
        pp.add("max_level", 0);
    }

    // Create the "input file"
    std::stringstream ss;
    ss << "2 // Number of levels" << '\n';
    ss << "2 // Number of boxes at this level" << '\n';
    ss << "-10.0 -75.0 0.0 15.0 -65.0 20.0" << '\n';
    ss << "-10.0 -55.0 0.0 15.0 -45.0 20.0" << '\n';
    ss << "2 // Number of boxes at this level" << '\n';
    ss << "-10.0  25.0 0.0 15.0  35.0 20.0" << '\n';
    ss << "-10.0  65.0 0.0 15.0  75.0 20.0" << '\n';

    {
        CaptureOutput io;
        create_mesh_instance<RefineMesh>();
        std::unique_ptr<amr_wind::CartBoxRefinement> box_refine(
            new amr_wind::CartBoxRefinement(sim()));
        box_refine->read_inputs(mesh(), ss);

        auto msg = io.stdout().str();
        EXPECT_GT(msg.size(), 0);
        auto found = msg.find("WARNING");
        EXPECT_NE(found, std::string::npos);
    }
}

/* Check that the implementation handles bounding box limits that extend beyond
 * the problem domain.
 */
TEST_F(NestRefineTest, bbox_limits)
{
    setup_refinement_inputs();

    // Create the "input file"
    std::stringstream ss;
    ss << "1 // Number of levels" << '\n';
    ss << "1 // Number of boxes at this level" << '\n';
    ss << "-60.0 -200.0 -10.0 35.0 200.0 60.0" << '\n';

    create_mesh_instance<RefineMesh>();
    std::unique_ptr<amr_wind::CartBoxRefinement> box_refine(
        new amr_wind::CartBoxRefinement(sim()));
    box_refine->read_inputs(mesh(), ss);

    auto targets = box_refine->boxarray_vec();
    EXPECT_EQ(targets.size(), 1U);
    EXPECT_EQ(targets[0].size(), 1U);

    auto domain = mesh().Geom(0).Domain();
    auto bx = targets[0][0];

    EXPECT_EQ(bx.smallEnd(), domain.smallEnd());
    auto big_end = domain.bigEnd();
    EXPECT_EQ(bx.bigEnd(), big_end.diagShift(1));
}

} // namespace amr_wind_tests
