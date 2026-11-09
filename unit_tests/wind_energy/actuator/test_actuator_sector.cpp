#include "ks_test_utils/MeshTest.H"

#include "src/wind_energy/actuator/Actuator.H"
#include "src/wind_energy/actuator/ActuatorContainer.H"
#include "src/wind_energy/actuator/ActuatorModel.H"
#include "src/wind_energy/actuator/ActParser.H"
#include "src/wind_energy/actuator/sector/ActuatorSector.H"
#include "src/wind_energy/actuator/sector/actuator_sector_ops.H"
#include "src/utilities/constants.H"
#include "src/utilities/tagging/ActuatorRefinement.H"

#include <fstream>
#include <limits>

#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {
namespace {

class ActuatorSectorTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{16, 16, 16}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 16);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", 1.0e-4_rt);
        }
        {
            amrex::ParmParse pp("geometry");
            pp.addarr(
                "prob_lo",
                amrex::Vector<amrex::Real>{-0.2_rt, -0.2_rt, -0.2_rt});
            pp.addarr(
                "prob_hi", amrex::Vector<amrex::Real>{0.2_rt, 0.2_rt, 0.2_rt});
        }
    }

    void initialize_domain()
    {
        initialize_mesh();
        sim().repo().declare_field("actuator_src_term", 3, 0);
        auto& vel = sim().repo().declare_field("velocity", 3, 3);
        auto& density = sim().repo().declare_field("density", 1, 3);
        vel.setVal(0.0_rt);
        density.setVal(1.0_rt);
        kynema_sgf::actuator::ActuatorContainer::ParticleType::NextID(1U);
    }

    void populate_actuator_inputs(
        const bool use_epsilon = false,
        const amrex::Real epsilon = 0.0_rt,
        const bool use_epsilon_chord = true,
        const amrex::Real epsilon_chord = 0.25_rt,
        const amrex::Real epsilon_min = 0.0_rt)
    {
        amrex::ParmParse pp_a("Actuator");
        pp_a.add("labels", (std::string) "R1");
        pp_a.add("type", (std::string) "ActuatorSector");

        amrex::ParmParse pp("Actuator.ActuatorSector");
        pp.add("rotor_diameter", 0.1_rt);
        pp.add("root_radius_fraction", 0.18_rt);
        pp.add("num_blades", 2);
        pp.add("omega", 2500.0_rt);
        if (use_epsilon) {
            pp.add("epsilon", epsilon);
        }
        if (use_epsilon_chord) {
            pp.add("epsilon_chord", epsilon_chord);
            pp.add("epsilon_min", epsilon_min);
        }
        pp.add("epsilon_dr", 1.5_rt);
        pp.add("epsilon_dl", 1.5_rt);
        pp.add("min_chord_dr", 2.0_rt);
        pp.add("airfoil_table", m_afname);
        pp.add("airfoil_type", (std::string) "openfast");
        pp.add("gaussian_type", (std::string) "table");
        pp.addarr("span_locs", amrex::Vector<amrex::Real>{0.0_rt, 1.0_rt});
        pp.addarr("chord", amrex::Vector<amrex::Real>{0.01_rt, 0.006_rt});
        pp.addarr("twist", amrex::Vector<amrex::Real>{12.0_rt, 4.0_rt});

        amrex::ParmParse pp_r("Actuator.R1");
        pp_r.addarr(
            "center", amrex::Vector<amrex::Real>{0.0_rt, 0.0_rt, 0.0_rt});
        pp_r.addarr(
            "translation_velocity",
            amrex::Vector<amrex::Real>{0.0_rt, 0.0_rt, 1.0_rt});
        pp_r.addarr(
            "rotor_normal", amrex::Vector<amrex::Real>{1.0_rt, 0.0_rt, 0.0_rt});
    }

    const std::string m_afname = "sector_airfoil.txt";
};

void write_airfoil_file(const std::string& fname)
{
    std::ofstream os(fname);
    os << "! OpenFAST-style test polar\n";
    os << "        5   NumAlf\n";
    os << "! Alpha Cl Cd Cm\n";
    os << "! deg   -  -  -\n";
    os << " -180.0  0.0 0.04 0.0\n";
    os << "  -10.0 -0.5 0.02 0.0\n";
    os << "    0.0  0.0 0.01 0.0\n";
    os << "   10.0  0.5 0.02 0.0\n";
    os << "  180.0  0.0 0.04 0.0\n";
}

void remove_file(const std::string& fname)
{
    std::ifstream f(fname);
    if (f.good()) {
        remove(fname.c_str());
    }
}

class ActPhysicsTest : public ::kynema_sgf::actuator::Actuator
{
public:
    explicit ActPhysicsTest(::kynema_sgf::CFDSim& sim)
        : ::kynema_sgf::actuator::Actuator(sim)
    {}

protected:
    void prepare_outputs() override {}
};

class ActuatorSectorRefinementTest : public ActuatorSectorTest
{
public:
    void populate_parameters() override
    {
        ActuatorSectorTest::populate_parameters();
        amrex::ParmParse pp("amr");
        pp.remove("max_level");
        pp.add("max_level", 1);
    }

    void create_mesh_instance() override
    {
        if (!m_mesh) {
            reset_prob_domain();
            m_mesh = std::make_unique<RefineMesh>();
        }
    }
};

} // namespace

TEST_F(ActuatorSectorTest, act_model_init)
{
    write_airfoil_file(m_afname);
    initialize_domain();
    populate_actuator_inputs();

    ::kynema_sgf::actuator::ActModel<
        ::kynema_sgf::actuator::ActuatorSector,
        ::kynema_sgf::actuator::ActSrcSector>
        sector(sim(), "R1", 0);
    {
        ::kynema_sgf::actuator::utils::ActParser pp(
            "Actuator.ActuatorSector", "Actuator.R1");
        sector.read_inputs(pp);
    }

    amrex::Vector<int> act_proc_count(amrex::ParallelDescriptor::NProcs(), 0);
    sector.determine_root_proc(act_proc_count);
    sector.init_actuator_source();

    const auto& meta = sector.meta();
    EXPECT_GT(meta.radius.size(), 0);
    EXPECT_EQ(meta.radius.size(), meta.dr.size());
    EXPECT_EQ(meta.gaussian_table_nintervals, 319);
    EXPECT_EQ(meta.gaussian_table.size(), meta.gaussian_table_nintervals + 1);
    EXPECT_EQ(sector.num_velocity_points(), 2 * meta.radius.size());

    const auto refinement = sector.refinement_geometries(0.5_rt);
    ASSERT_EQ(refinement.size(), 1U);
    EXPECT_EQ(refinement[0].label, "R1");
    constexpr amrex::Real tol = kynema_sgf::constants::TIGHT_TOL;
    EXPECT_NEAR(refinement[0].center.x(), 0.0_rt, tol);
    EXPECT_NEAR(refinement[0].center.y(), 0.0_rt, tol);
    EXPECT_NEAR(refinement[0].center.z(), 0.5_rt, tol);
    EXPECT_NEAR(refinement[0].normal.x(), 1.0_rt, tol);
    EXPECT_NEAR(refinement[0].normal.y(), 0.0_rt, tol);
    EXPECT_NEAR(refinement[0].normal.z(), 0.0_rt, tol);
    EXPECT_NEAR(refinement[0].rotor_radius, 0.05_rt, tol);
    EXPECT_NEAR(refinement[0].epsilon_max, 0.0025_rt, tol);

    remove_file(m_afname);
}

TEST_F(ActuatorSectorTest, actuator_lifecycle)
{
    write_airfoil_file(m_afname);
    initialize_domain();
    populate_actuator_inputs();

    ActPhysicsTest act(sim());
    act.pre_init_actions();
    act.post_init_actions();
    act.pre_advance_work();

    remove_file(m_afname);
}

TEST_F(ActuatorSectorRefinementTest, asymmetric_actuator_refinement)
{
    write_airfoil_file(m_afname);
    initialize_domain();
    populate_actuator_inputs();

    auto& physics = sim().physics_manager().create("Actuator", sim());
    auto& actuator = dynamic_cast<::kynema_sgf::actuator::Actuator&>(physics);
    actuator.pre_init_actions();

    amrex::ParmParse pp("tagging.rotor_tracking");
    pp.addarr("actuator_labels", amrex::Vector<std::string>{"R1"});
    pp.add("min_level", 0);
    pp.add("max_level", 0);
    pp.add("radial_padding_epsilon", 0.0_rt);
    pp.add("forward_padding_diameter", 0.8_rt);
    pp.add("backward_padding_diameter", 0.3_rt);

    auto refiner =
        ::kynema_sgf::RefinementCriteria::create("ActuatorRefinement", sim());
    refiner->initialize("tagging.rotor_tracking");

    const auto& field = sim().repo().get_field("velocity")(0);
    amrex::TagBoxArray tags(field.boxArray(), field.DistributionMap());
    tags.setVal(amrex::TagBox::CLEAR);
    (*refiner)(0, tags, 0.0_rt, 0);
    amrex::Gpu::streamSynchronize();

    amrex::Gpu::PinnedVector<amrex::IntVect> tagged_cells;
    tags.collate(tagged_cells);
    ASSERT_EQ(tagged_cells.size(), 48U);

    const auto& geom = sim().mesh().Geom(0);
    const auto& problo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();
    amrex::Real xmin = std::numeric_limits<amrex::Real>::max();
    amrex::Real xmax = std::numeric_limits<amrex::Real>::lowest();
    for (const auto& cell : tagged_cells) {
        const amrex::Real x = problo[0] + ((cell[0] + 0.5_rt) * dx[0]);
        xmin = amrex::min(xmin, x);
        xmax = amrex::max(xmax, x);
    }
    constexpr amrex::Real tol = kynema_sgf::constants::TIGHT_TOL;
    EXPECT_NEAR(xmin, -0.0125_rt, tol);
    EXPECT_NEAR(xmax, 0.0625_rt, tol);

    remove_file(m_afname);
}

TEST_F(ActuatorSectorRefinementTest, logical_operators_use_cylinder_union)
{
    write_airfoil_file(m_afname);
    initialize_domain();
    populate_actuator_inputs();

    amrex::ParmParse pp_a("Actuator");
    pp_a.addarr("labels", amrex::Vector<std::string>{"R1", "R2"});
    amrex::ParmParse pp_r("Actuator.R2");
    pp_r.addarr("center", amrex::Vector<amrex::Real>{0.05_rt, 0.0_rt, 0.0_rt});
    pp_r.addarr(
        "rotor_normal", amrex::Vector<amrex::Real>{1.0_rt, 0.0_rt, 0.0_rt});

    auto& physics = sim().physics_manager().create("Actuator", sim());
    auto& actuator = dynamic_cast<::kynema_sgf::actuator::Actuator&>(physics);
    actuator.pre_init_actions();

    amrex::ParmParse pp("tagging.rotor_tracking");
    pp.addarr("actuator_labels", amrex::Vector<std::string>{"R1", "R2"});
    pp.add("radial_padding_epsilon", 0.0_rt);
    pp.add("forward_padding_diameter", 0.8_rt);
    pp.add("backward_padding_diameter", 0.3_rt);

    const auto& field = sim().repo().get_field("velocity")(0);
    amrex::TagBoxArray tags(field.boxArray(), field.DistributionMap());
    const auto& geom = sim().mesh().Geom(0);
    const auto& problo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();
    const auto total_cells = static_cast<std::size_t>(geom.Domain().numPts());
    constexpr std::size_t union_cells = 72;

    for (const std::string op : {"or", "and", "and_not", "or_not"}) {
        amrex::ParmParse pp_tagging("tagging");
        pp_tagging.addarr(
            "labels", amrex::Vector<std::string>{"rotor_tracking"});
        pp.add("type", std::string("ActuatorRefinement"));
        pp.add("operator", op);
        ::kynema_sgf::RefineCriteriaManager manager(sim());
        manager.initialize();
        for (const bool previous_tag : {false, true}) {
            SCOPED_TRACE(
                op + (previous_tag ? " initially set" : " initially clear"));
            tags.setVal(
                previous_tag ? amrex::TagBox::SET : amrex::TagBox::CLEAR);
            manager.tag_cells(1, tags, 0.0_rt, 0);
            amrex::Gpu::PinnedVector<amrex::IntVect> inactive_cells;
            tags.collate(inactive_cells);
            EXPECT_EQ(inactive_cells.size(), previous_tag ? total_cells : 0);
            tags.setVal(
                previous_tag ? amrex::TagBox::SET : amrex::TagBox::CLEAR);
            manager.tag_cells(0, tags, 0.0_rt, 0);
            amrex::Gpu::streamSynchronize();

            amrex::Gpu::PinnedVector<amrex::IntVect> tagged_cells;
            tags.collate(tagged_cells);
            std::size_t expected_count = 0;
            if (op == "or") {
                expected_count = previous_tag ? total_cells : union_cells;
            } else if (op == "and") {
                expected_count = previous_tag ? union_cells : 0;
            } else if (op == "and_not") {
                expected_count = previous_tag ? total_cells - union_cells : 0;
            } else {
                expected_count =
                    previous_tag ? total_cells : total_cells - union_cells;
            }
            EXPECT_EQ(tagged_cells.size(), expected_count);

            for (const auto& cell : tagged_cells) {
                const amrex::Real xpos =
                    problo[0] + ((cell[0] + 0.5_rt) * dx[0]);
                const amrex::Real ypos =
                    problo[1] + ((cell[1] + 0.5_rt) * dx[1]);
                const amrex::Real zpos =
                    problo[2] + ((cell[2] + 0.5_rt) * dx[2]);
                const bool in_union =
                    ((xpos >= -0.03_rt && xpos <= 0.08_rt) ||
                     (xpos >= 0.02_rt && xpos <= 0.13_rt)) &&
                    ((ypos * ypos) + (zpos * zpos)) <= 0.05_rt * 0.05_rt;
                bool expected_tag = false;
                if (op == "or") {
                    expected_tag = previous_tag || in_union;
                } else if (op == "and") {
                    expected_tag = previous_tag && in_union;
                } else if (op == "and_not") {
                    expected_tag = previous_tag && !in_union;
                } else {
                    expected_tag = previous_tag || !in_union;
                }
                EXPECT_TRUE(expected_tag);
            }
        }
    }

    remove_file(m_afname);
}

TEST_F(ActuatorSectorTest, fixed_epsilon_input)
{
    write_airfoil_file(m_afname);
    initialize_domain();
    populate_actuator_inputs(true, 0.004_rt, false);

    ::kynema_sgf::actuator::ActModel<
        ::kynema_sgf::actuator::ActuatorSector,
        ::kynema_sgf::actuator::ActSrcSector>
        sector(sim(), "R1", 0);
    {
        ::kynema_sgf::actuator::utils::ActParser pp(
            "Actuator.ActuatorSector", "Actuator.R1");
        sector.read_inputs(pp);
    }
    sector.init_actuator_source();

    const auto& meta = sector.meta();
    for (const auto eps : meta.epsilon_profile) {
        EXPECT_NEAR(eps, 0.004_rt, 1.0e-14_rt);
    }

    remove_file(m_afname);
}

TEST_F(ActuatorSectorTest, epsilon_chord_uses_minimum)
{
    write_airfoil_file(m_afname);
    initialize_domain();
    populate_actuator_inputs(false, 0.0_rt, true, 0.25_rt, 0.004_rt);

    ::kynema_sgf::actuator::ActModel<
        ::kynema_sgf::actuator::ActuatorSector,
        ::kynema_sgf::actuator::ActSrcSector>
        sector(sim(), "R1", 0);
    {
        ::kynema_sgf::actuator::utils::ActParser pp(
            "Actuator.ActuatorSector", "Actuator.R1");
        sector.read_inputs(pp);
    }
    sector.init_actuator_source();

    const auto& meta = sector.meta();
    for (int i = 0; i < static_cast<int>(meta.epsilon_profile.size()); ++i) {
        const amrex::Real expected =
            amrex::max(0.004_rt, 0.25_rt * meta.chord[i]);
        EXPECT_NEAR(meta.epsilon_profile[i], expected, 1.0e-14_rt);
    }

    remove_file(m_afname);
}

} // namespace kynema_sgf_tests
