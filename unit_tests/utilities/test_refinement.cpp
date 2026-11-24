#include <sstream>

#include "ks_test_utils/AmrexTest.H"
#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/OutputCapture.H"

#include "AMReX_Box.H"
#include "AMReX_BoxArray.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_BoxList.H"
#include "AMReX_Geometry.H"
#include "AMReX_RealBox.H"
#include "AMReX_Vector.H"

#include "src/utilities/tagging/CartBoxRefinement.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

template <typename Predicate>
void expect_tag_predicate(
    amrex::TagBoxArray& tags,
    const amrex::Box& domain,
    const Predicate& predicate)
{
    amrex::Gpu::PinnedVector<amrex::IntVect> cells;
    tags.collate(cells);
    std::size_t expected_count = 0;
    for (amrex::BoxIterator cell(domain); cell.ok(); ++cell) {
        if (predicate(cell())) {
            ++expected_count;
        }
    }
    EXPECT_EQ(cells.size(), expected_count);
    for (const auto& cell : cells) {
        EXPECT_TRUE(predicate(cell));
    }
}

void tag_parity(amrex::TagBoxArray& tags, const int parity)
{
    for (amrex::MFIter mfi(tags); mfi.isValid(); ++mfi) {
        const auto& tag = tags.array(mfi);
        amrex::ParallelFor(
            mfi.validbox(), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (i % 2 == parity) {
                    tag(i, j, k) = amrex::TagBox::SET;
                }
            });
    }
}

class DummyOperatorRefinement
    : public kynema_sgf::RefinementCriteria::Register<DummyOperatorRefinement>
{
public:
    static std::string identifier() { return "DummyOperatorRefinement"; }

    explicit DummyOperatorRefinement(kynema_sgf::CFDSim& /*sim*/) {}

    void initialize(const std::string& key) override
    {
        amrex::ParmParse pp(key);
        pp.query("parity", m_parity);
        pp.query("max_level", m_max_level);
    }

    [[nodiscard]] bool active(int level) const override
    {
        return level <= m_max_level;
    }

    void operator()(
        int /*level*/,
        amrex::TagBoxArray& tags,
        amrex::Real /*time*/,
        int /*ngrow*/) override
    {
        tag_parity(tags, m_parity);
    }

private:
    int m_parity{0};
    int m_max_level{32};
};

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
    using kynema_sgf::tagging::string_to_operator;
    using kynema_sgf::tagging::TaggingOperator;

    EXPECT_EQ(string_to_operator("and"), TaggingOperator::AND);
    EXPECT_EQ(string_to_operator("or"), TaggingOperator::OR);
    EXPECT_EQ(string_to_operator("and_not"), TaggingOperator::AND_NOT);
    EXPECT_EQ(string_to_operator("or_not"), TaggingOperator::OR_NOT);
}

TEST(RefinementTaggingLogic, tag_val_truth_table)
{
    using kynema_sgf::tagging::tag_val;
    using kynema_sgf::tagging::TaggingOperator;

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

TEST_F(NestRefineTest, manager_applies_operator_per_label)
{
    setup_refinement_inputs();
    create_mesh_instance<RefineMesh>();
    initialize_mesh();
    amrex::TagBoxArray tags(mesh().boxArray(0), mesh().DistributionMap(0));
    amrex::ParmParse pp_tagging("tagging");
    pp_tagging.addarr("labels", amrex::Vector<std::string>{"criterion"});
    amrex::ParmParse pp("tagging.criterion");
    pp.add("type", std::string("DummyOperatorRefinement"));

    for (const std::string input_op : {"and", "OR", "and_not", "Or_Not", ""}) {
        pp.remove("operator");
        if (!input_op.empty()) {
            pp.add("operator", input_op);
        }
        const auto op = input_op.empty() ? "or" : amrex::toLower(input_op);
        kynema_sgf::RefineCriteriaManager manager(sim());
        manager.initialize();
        for (const bool previous : {false, true}) {
            SCOPED_TRACE(input_op + (previous ? " set" : " clear"));
            tags.setVal(previous ? amrex::TagBox::SET : amrex::TagBox::CLEAR);
            manager.tag_cells(0, tags, 0.0_rt, 0);
            expect_tag_predicate(
                tags, mesh().Geom(0).Domain(), [=](const amrex::IntVect& cell) {
                    const bool current = cell[0] % 2 == 0;
                    if (op == "and") {
                        return previous && current;
                    }
                    if (op == "and_not") {
                        return previous && !current;
                    }
                    if (op == "or_not") {
                        return previous || !current;
                    }
                    return previous || current;
                });
        }
    }
}

TEST_F(NestRefineTest, manager_preserves_order_and_clears_scratch)
{
    setup_refinement_inputs();
    create_mesh_instance<RefineMesh>();
    initialize_mesh();
    amrex::TagBoxArray tags(mesh().boxArray(0), mesh().DistributionMap(0));
    amrex::ParmParse pp_even("tagging.even");
    pp_even.add("type", std::string("DummyOperatorRefinement"));
    pp_even.add("operator", std::string("or"));
    pp_even.add("parity", 0);
    amrex::ParmParse pp_odd("tagging.odd");
    pp_odd.add("type", std::string("DummyOperatorRefinement"));
    pp_odd.add("operator", std::string("and"));
    pp_odd.add("parity", 1);

    for (const bool reverse : {false, true}) {
        SCOPED_TRACE(reverse ? "odd then even" : "even then odd");
        amrex::ParmParse pp("tagging");
        pp.addarr(
            "labels", reverse ? amrex::Vector<std::string>{"odd", "even"}
                              : amrex::Vector<std::string>{"even", "odd"});
        kynema_sgf::RefineCriteriaManager manager(sim());
        manager.initialize();
        tags.setVal(amrex::TagBox::CLEAR);
        manager.tag_cells(0, tags, 0.0_rt, 0);
        expect_tag_predicate(
            tags, mesh().Geom(0).Domain(), [=](const amrex::IntVect& cell) {
                return reverse && cell[0] % 2 == 0;
            });
    }
}

TEST_F(NestRefineTest, manager_distinguishes_inactive_and_empty)
{
    setup_refinement_inputs();
    create_mesh_instance<RefineMesh>();
    initialize_mesh();
    amrex::TagBoxArray tags(mesh().boxArray(0), mesh().DistributionMap(0));
    amrex::ParmParse pp_tagging("tagging");
    pp_tagging.addarr("labels", amrex::Vector<std::string>{"empty"});
    amrex::ParmParse pp("tagging.empty");
    pp.add("type", std::string("DummyOperatorRefinement"));
    pp.add("parity", -1);

    for (const std::string op : {"and", "or_not"}) {
        pp.add("operator", op);
        for (const bool active : {false, true}) {
            SCOPED_TRACE(op + (active ? " active" : " inactive"));
            pp.add("max_level", active ? 0 : -1);
            kynema_sgf::RefineCriteriaManager manager(sim());
            manager.initialize();
            const bool previous = op == "and";
            tags.setVal(previous ? amrex::TagBox::SET : amrex::TagBox::CLEAR);
            manager.tag_cells(0, tags, 0.0_rt, 0);
            expect_tag_predicate(
                tags, mesh().Geom(0).Domain(), [=](const amrex::IntVect&) {
                    return active ? !previous : previous;
                });
        }
    }

    pp_tagging.addarr("labels", amrex::Vector<std::string>{});
    kynema_sgf::RefineCriteriaManager manager(sim());
    manager.initialize();
    tags.setVal(amrex::TagBox::SET);
    manager.tag_cells(0, tags, 0.0_rt, 0);
    expect_tag_predicate(
        tags, mesh().Geom(0).Domain(),
        [](const amrex::IntVect&) { return true; });
}

TEST_F(NestRefineTest, manager_combines_geometry_union_once)
{
    setup_refinement_inputs();
    create_mesh_instance<RefineMesh>();
    initialize_mesh();
    sim().repo().declare_field("tagging_test", 1, 0);
    amrex::TagBoxArray tags(mesh().boxArray(0), mesh().DistributionMap(0));
    amrex::ParmParse pp_tagging("tagging");
    pp_tagging.addarr("labels", amrex::Vector<std::string>{"geometry"});
    amrex::ParmParse pp("tagging.geometry");
    pp.add("type", std::string("GeometryRefinement"));
    pp.addarr("shapes", amrex::Vector<std::string>{"left", "right"});
    pp.add("max_level", 0);
    for (const std::string shape : {"left", "right"}) {
        amrex::ParmParse pp_shape("tagging.geometry." + shape);
        pp_shape.add("type", std::string("cylinder"));
        const bool left = shape == "left";
        pp_shape.addarr(
            "start", amrex::Vector<amrex::Real>{
                         left ? -15.0_rt : -5.0_rt, 0.0_rt, 15.0_rt});
        pp_shape.addarr(
            "end", amrex::Vector<amrex::Real>{
                       left ? 5.0_rt : 15.0_rt, 0.0_rt, 15.0_rt});
        pp_shape.add("outer_radius", 10.0_rt);
    }
    const auto& geom = mesh().Geom(0);
    const auto problo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();
    for (const std::string op : {"and", "and_not", "or_not"}) {
        pp.add("operator", op);
        kynema_sgf::RefineCriteriaManager manager(sim());
        manager.initialize();
        const bool previous = op != "or_not";
        tags.setVal(previous ? amrex::TagBox::SET : amrex::TagBox::CLEAR);
        manager.tag_cells(0, tags, 0.0_rt, 0);
        expect_tag_predicate(
            tags, geom.Domain(), [=](const amrex::IntVect& cell) {
                const auto xpos = problo[0] + ((cell[0] + 0.5_rt) * dx[0]);
                const auto ypos = problo[1] + ((cell[1] + 0.5_rt) * dx[1]);
                const auto zpos =
                    problo[2] + ((cell[2] + 0.5_rt) * dx[2]) - 15.0_rt;
                const bool inside = xpos >= -15.0_rt && xpos <= 15.0_rt &&
                                    (ypos * ypos) + (zpos * zpos) <= 100.0_rt;
                return op == "and" ? inside : !inside;
            });
    }
}

TEST_F(NestRefineTest, manager_combines_field_and_gradient_union_once)
{
    setup_refinement_inputs();
    create_mesh_instance<RefineMesh>();
    initialize_mesh();
    auto& field = sim().repo().declare_int_field("indicator", 1, 1);
    field.setVal(1);
    amrex::TagBoxArray tags(mesh().boxArray(0), mesh().DistributionMap(0));
    amrex::ParmParse pp_tagging("tagging");
    pp_tagging.addarr("labels", amrex::Vector<std::string>{"field"});
    amrex::ParmParse pp("tagging.field");
    pp.add("type", std::string("FieldRefinement"));
    pp.add("field_name", std::string("indicator"));
    pp.addarr("field_error", amrex::Vector<amrex::Real>{0.5_rt});
    pp.addarr("grad_error", amrex::Vector<amrex::Real>{1.0_rt});
    for (const std::string op : {"and", "or_not"}) {
        pp.add("operator", op);
        kynema_sgf::RefineCriteriaManager manager(sim());
        manager.initialize();
        const bool previous = op == "and";
        tags.setVal(previous ? amrex::TagBox::SET : amrex::TagBox::CLEAR);
        manager.tag_cells(0, tags, 0.0_rt, 0);
        expect_tag_predicate(
            tags, mesh().Geom(0).Domain(),
            [=](const amrex::IntVect&) { return previous; });
        manager.tag_cells(1, tags, 0.0_rt, 0);
        expect_tag_predicate(
            tags, mesh().Geom(0).Domain(),
            [=](const amrex::IntVect&) { return previous; });
    }
}

} // namespace kynema_sgf_tests
