#include <sstream>

#include "ks_test_utils/AmrexTest.H"
#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/OutputCapture.H"

#include "AMReX_Box.H"
#include "AMReX_BoxArray.H"
#include "AMReX_BoxList.H"
#include "AMReX_Geometry.H"
#include "AMReX_RealBox.H"
#include "AMReX_Vector.H"

#include "src/utilities/tagging/CartBoxRefinement.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

class DummyOperatorRefinement
    : public kynema_sgf::RefinementCriteria::Register<DummyOperatorRefinement>
{
public:
    static std::string identifier() { return "DummyOperatorRefinement"; }

    explicit DummyOperatorRefinement(kynema_sgf::CFDSim& /*sim*/) {}

    static void reset_seen_ops() { s_seen_ops.clear(); }

    static const amrex::Vector<kynema_sgf::tagging::TaggingOperator>& seen_ops()
    {
        return s_seen_ops;
    }

    void initialize(const std::string& /*key*/) override
    {
        s_seen_ops.push_back(tag_operator());
    }

    void operator()(
        int /*level*/,
        amrex::TagBoxArray& /*tags*/,
        amrex::Real /*time*/,
        int /*ngrow*/) override
    {}

private:
    static amrex::Vector<kynema_sgf::tagging::TaggingOperator> s_seen_ops;
};

amrex::Vector<kynema_sgf::tagging::TaggingOperator>
    DummyOperatorRefinement::s_seen_ops;

} // namespace

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
    ref_vec.emplace_back(
        std::make_unique<kynema_sgf::CartBoxRefinement>(sim()));
    auto* box_refine =
        dynamic_cast<kynema_sgf::CartBoxRefinement*>(ref_vec[0].get());
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
        std::unique_ptr<kynema_sgf::CartBoxRefinement> box_refine(
            new kynema_sgf::CartBoxRefinement(sim()));
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
    std::unique_ptr<kynema_sgf::CartBoxRefinement> box_refine(
        new kynema_sgf::CartBoxRefinement(sim()));
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

TEST(RefinementTaggingLogic, string_to_operator)
{
    using kynema_sgf::tagging::TaggingOperator;
    using kynema_sgf::tagging::string_to_operator;

    EXPECT_EQ(string_to_operator("and"), TaggingOperator::AND);
    EXPECT_EQ(string_to_operator("OR"), TaggingOperator::OR);
    EXPECT_EQ(string_to_operator("And_Not"), TaggingOperator::AND_NOT);
    EXPECT_EQ(string_to_operator("or_not"), TaggingOperator::OR_NOT);
}

TEST(RefinementTaggingLogic, tag_val_truth_table)
{
    using kynema_sgf::tagging::TaggingOperator;
    using kynema_sgf::tagging::tag_val;

    auto is_set = [](amrex::TagBox::TagVal val) {
        return val == amrex::TagBox::SET;
    };

    // and
    EXPECT_FALSE(is_set(tag_val(false, false, TaggingOperator::AND)));
    EXPECT_FALSE(is_set(tag_val(false, true, TaggingOperator::AND)));
    EXPECT_FALSE(is_set(tag_val(true, false, TaggingOperator::AND)));
    EXPECT_TRUE(is_set(tag_val(true, true, TaggingOperator::AND)));

    // or
    EXPECT_FALSE(is_set(tag_val(false, false, TaggingOperator::OR)));
    EXPECT_TRUE(is_set(tag_val(false, true, TaggingOperator::OR)));
    EXPECT_TRUE(is_set(tag_val(true, false, TaggingOperator::OR)));
    EXPECT_TRUE(is_set(tag_val(true, true, TaggingOperator::OR)));

    // and_not
    EXPECT_FALSE(is_set(tag_val(false, false, TaggingOperator::AND_NOT)));
    EXPECT_FALSE(is_set(tag_val(false, true, TaggingOperator::AND_NOT)));
    EXPECT_TRUE(is_set(tag_val(true, false, TaggingOperator::AND_NOT)));
    EXPECT_FALSE(is_set(tag_val(true, true, TaggingOperator::AND_NOT)));

    // or_not
    EXPECT_TRUE(is_set(tag_val(false, false, TaggingOperator::OR_NOT)));
    EXPECT_FALSE(is_set(tag_val(false, true, TaggingOperator::OR_NOT)));
    EXPECT_TRUE(is_set(tag_val(true, false, TaggingOperator::OR_NOT)));
    EXPECT_TRUE(is_set(tag_val(true, true, TaggingOperator::OR_NOT)));
}

TEST_F(NestRefineTest, manager_parses_operator_per_label)
{
    setup_refinement_inputs();
    create_mesh_instance<RefineMesh>();

    DummyOperatorRefinement::reset_seen_ops();

    {
        amrex::ParmParse pp("tagging");
        amrex::Vector<std::string> labels{{"t_and", "t_or", "t_and_not", "t_or_not"}};
        pp.addarr("labels", labels);
    }
    {
        amrex::ParmParse pp("tagging.t_and");
        pp.add("type", (std::string) "DummyOperatorRefinement");
        pp.add("operator", (std::string) "and");
    }
    {
        amrex::ParmParse pp("tagging.t_or");
        pp.add("type", (std::string) "DummyOperatorRefinement");
        pp.add("operator", (std::string) "OR");
    }
    {
        amrex::ParmParse pp("tagging.t_and_not");
        pp.add("type", (std::string) "DummyOperatorRefinement");
        pp.add("operator", (std::string) "and_not");
    }
    {
        amrex::ParmParse pp("tagging.t_or_not");
        pp.add("type", (std::string) "DummyOperatorRefinement");
        pp.add("operator", (std::string) "Or_Not");
    }

    kynema_sgf::RefineCriteriaManager manager(sim());
    manager.initialize();

    const auto& ops = DummyOperatorRefinement::seen_ops();
    ASSERT_EQ(ops.size(), 4U);
    EXPECT_EQ(ops[0], kynema_sgf::tagging::TaggingOperator::AND);
    EXPECT_EQ(ops[1], kynema_sgf::tagging::TaggingOperator::OR);
    EXPECT_EQ(ops[2], kynema_sgf::tagging::TaggingOperator::AND_NOT);
    EXPECT_EQ(ops[3], kynema_sgf::tagging::TaggingOperator::OR_NOT);
}

} // namespace kynema_sgf_tests
