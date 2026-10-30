#include "ks_test_utils/MeshTest.H"
#include "src/boundary_conditions/BCInterface.H"
#include "src/boundary_conditions/field_boundary_fill/BoundaryPlane.H"
#include "src/boundary_conditions/field_boundary_fill/ModulatedPowerLaw.H"
#include "src/equation_systems/vof/vof_bciface.H"

namespace kynema_sgf_tests {

namespace {
const amrex::Orientation xlo(0, amrex::Orientation::low);
const amrex::Orientation xhi(0, amrex::Orientation::high);
} // namespace

class BCInputTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        amrex::ParmParse("geometry")
            .addarr("is_periodic", amrex::Vector<int>{0, 0, 0});
        for (const auto& face : kynema_sgf::bcnames) {
            amrex::ParmParse(face).add("type", std::string("zero_gradient"));
        }
        amrex::ParmParse("xlo").add("type", std::string("mass_inflow"));
    }

    kynema_sgf::Field& scalar(const std::string& name = "temperature")
    {
        auto& field = sim().repo().declare_field(name, 1, 1);
        kynema_sgf::BCScalar bc(field);
        bc();
        return field;
    }

    void validate()
    {
        kynema_sgf::BCIface::validate_inflow_values(
            sim().repo(), sim().field_boundary_manager());
    }

    kynema_sgf::BoundaryPlane& plane(
        const amrex::Vector<std::string>& fields,
        const std::string& fluid_phase = "both",
        const bool input = true,
        const bool legacy = false)
    {
        amrex::ParmParse pp(legacy ? "ABL" : "BoundaryPlane");
        const std::string prefix = legacy ? "bndry_" : "";
        pp.add(prefix + "io_mode", input ? 1 : 0);
        pp.add(prefix + "file", std::string("unused"));
        pp.addarr(prefix + "var_names", fields);
        pp.add(prefix + "fluid_phase", fluid_phase);
        // Header setup in this format needs no on-disk input or ERF reader.
        pp.add(prefix + "output_format", std::string("erf-multiblock"));
        if (!input) {
            pp.addarr(prefix + "planes", amrex::Vector<std::string>{"xlo"});
        }
        auto& boundary = dynamic_cast<kynema_sgf::BoundaryPlane&>(
            sim().field_boundary_manager().create("BoundaryPlane", sim()));
        if (input) {
            boundary.initialize_data();
            boundary.read_header();
        }
        return boundary;
    }
};

TEST_F(BCInputTest, missing_constants_are_rejected_for_velocity_and_scalars)
{
    initialize_mesh();
    auto& velocity = sim().repo().declare_field("velocity", 3, 1);
    kynema_sgf::BCVelocity bc(velocity);
    bc();
    EXPECT_THROW(validate(), amrex::RuntimeError);

    amrex::ParmParse("xlo").addarr(
        "velocity", amrex::Vector<amrex::Real>{1.0, 2.0, 3.0});
    EXPECT_NO_THROW(validate());
    EXPECT_NO_THROW(scalar());
    EXPECT_THROW(validate(), amrex::RuntimeError);

    amrex::ParmParse("xlo").addarr(
        "temperature", amrex::Vector<amrex::Real>{300.0});
    EXPECT_NO_THROW(validate());
}

TEST_F(BCInputTest, constants_are_preserved)
{
    initialize_mesh();
    amrex::ParmParse("xlo").addarr(
        "velocity", amrex::Vector<amrex::Real>{1.0, 2.0, 3.0});
    amrex::ParmParse("xlo").addarr(
        "temperature", amrex::Vector<amrex::Real>{300.0});
    auto& velocity = sim().repo().declare_field("velocity", 3, 1);
    kynema_sgf::BCVelocity bc(velocity);
    bc();
    auto& temperature = scalar();
    EXPECT_NO_THROW(validate());
    EXPECT_EQ(
        velocity.bc_values()[xlo], (amrex::Vector<amrex::Real>{1.0, 2.0, 3.0}));
    EXPECT_EQ(temperature.bc_values()[xlo][0], 300.0);
}

TEST_F(BCInputTest, profiles_supply_inflow_and_inflow_outflow_values)
{
    initialize_mesh();
    amrex::ParmParse("xlo").add(
        "temperature.inflow_type", std::string("LinearProfile"));
    amrex::ParmParse("xhi").add("type", std::string("mass_inflow_outflow"));
    amrex::ParmParse("xhi").add(
        "temperature.inflow_outflow_type", std::string("LinearProfile"));
    scalar();
    EXPECT_NO_THROW(validate());
}

TEST_F(BCInputTest, inflow_outflow_constants_are_required)
{
    initialize_mesh();
    amrex::ParmParse("xlo").add("type", std::string("mass_inflow_outflow"));
    auto& velocity = sim().repo().declare_field("velocity", 3, 1);
    kynema_sgf::BCVelocity bc(velocity);
    bc();
    EXPECT_THROW(validate(), amrex::RuntimeError);
    amrex::ParmParse("xlo").addarr(
        "velocity", amrex::Vector<amrex::Real>{0.0, 0.0, 0.0});
    scalar();
    EXPECT_THROW(validate(), amrex::RuntimeError);
    plane({"temperature"});
    EXPECT_NO_THROW(validate());
}

TEST_F(BCInputTest, scalar_defaults_do_not_replace_required_inputs)
{
    initialize_mesh();
    auto& density = sim().repo().declare_field("density", 1, 1);
    kynema_sgf::BCScalar bc(density);
    bc(1.0);
    EXPECT_THROW(validate(), amrex::RuntimeError);
}

TEST_F(BCInputTest, non_inflow_defaults_are_preserved)
{
    initialize_mesh();
    amrex::ParmParse("xlo").add("type", std::string("no_slip_wall"));
    auto& velocity = sim().repo().declare_field("velocity", 3, 1);
    kynema_sgf::BCVelocity bc(velocity);
    bc();
    auto& temperature = scalar();
    EXPECT_NO_THROW(validate());
    EXPECT_EQ(
        velocity.bc_values()[xlo], (amrex::Vector<amrex::Real>{0.0, 0.0, 0.0}));
    EXPECT_EQ(temperature.bc_values()[xlo][0], 0.0);
}

TEST_F(BCInputTest, full_plane_supplies_only_selected_fields)
{
    initialize_mesh();
    auto& temperature = scalar();
    auto& boundary = plane({"temperature"});
    EXPECT_TRUE(boundary.provides_inflow_values(temperature, xlo));
    EXPECT_NO_THROW(validate());
    scalar("density");
    EXPECT_THROW(validate(), amrex::RuntimeError);
}

TEST_F(BCInputTest, plane_must_supply_the_face)
{
    initialize_mesh();
    auto& temperature = scalar();
    auto& boundary = plane({"temperature"});
    EXPECT_FALSE(boundary.provides_inflow_values(temperature, xhi));
    amrex::ParmParse("xhi").add("type", std::string("mass_inflow_outflow"));
    kynema_sgf::BCScalar bc(temperature);
    bc();
    EXPECT_THROW(validate(), amrex::RuntimeError);
}

TEST_F(BCInputTest, legacy_plane_inputs_supply_values)
{
    initialize_mesh();
    scalar();
    plane({"temperature"}, "both", true, true);
    EXPECT_NO_THROW(validate());
}

TEST_F(BCInputTest, output_planes_do_not_supply_values)
{
    initialize_mesh();
    auto& temperature = scalar();
    auto& boundary = plane({"temperature"}, "both", false);
    EXPECT_FALSE(boundary.provides_inflow_values(temperature, xlo));
    EXPECT_THROW(validate(), amrex::RuntimeError);
}

TEST_F(BCInputTest, phase_limited_planes_require_a_fallback)
{
    initialize_mesh();
    auto& temperature = scalar();
    auto& boundary = plane({"temperature"}, "gas");
    EXPECT_FALSE(boundary.provides_inflow_values(temperature, xlo));
    EXPECT_THROW(validate(), amrex::RuntimeError);
    amrex::ParmParse("xlo").addarr(
        "temperature", amrex::Vector<amrex::Real>{300.0});
    testing::internal::CaptureStdout();
    EXPECT_NO_THROW(validate());
    EXPECT_EQ(
        testing::internal::GetCapturedStdout().find("Warning:"),
        std::string::npos);
}

TEST_F(BCInputTest, phase_limited_planes_can_use_a_profile_fallback)
{
    initialize_mesh();
    amrex::ParmParse("xlo").add(
        "temperature.inflow_type", std::string("LinearProfile"));
    scalar();
    plane({"temperature"}, "liquid");
    EXPECT_NO_THROW(validate());
}

TEST_F(BCInputTest, phase_limited_velocity_planes_require_a_fallback)
{
    initialize_mesh();
    auto& velocity = sim().repo().declare_field("velocity", 3, 1);
    kynema_sgf::BCVelocity bc(velocity);
    bc();
    plane({"velocity"}, "liquid");
    EXPECT_THROW(validate(), amrex::RuntimeError);
}

TEST_F(BCInputTest, phase_restriction_does_not_apply_to_density_or_vof)
{
    initialize_mesh();
    auto& density = sim().repo().declare_field("density", 1, 1);
    auto& vof = sim().repo().declare_field("vof", 1, 1);
    kynema_sgf::BCDensityMod bc_density(density);
    kynema_sgf::BCVOF bc_vof(vof);
    bc_density();
    bc_vof();
    EXPECT_THROW(validate(), amrex::RuntimeError);
    auto& boundary = plane({"density", "vof"}, "liquid");
    EXPECT_TRUE(boundary.provides_inflow_values(density, xlo));
    EXPECT_TRUE(boundary.provides_inflow_values(vof, xlo));
    EXPECT_NO_THROW(validate());
}

TEST_F(BCInputTest, redundant_constants_produce_a_warning)
{
    initialize_mesh();
    amrex::ParmParse("xlo").addarr(
        "temperature", amrex::Vector<amrex::Real>{300.0});
    scalar();
    plane({"temperature"});
    testing::internal::CaptureStdout();
    EXPECT_NO_THROW(validate());
    const auto output = testing::internal::GetCapturedStdout();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        EXPECT_NE(
            output.find("Warning: xlo.temperature is overridden"),
            std::string::npos);
    }
}

TEST_F(BCInputTest, profile_constants_produce_a_warning_once_per_face)
{
    initialize_mesh();
    amrex::ParmParse("xlo").add(
        "temperature.inflow_type", std::string("LinearProfile"));
    amrex::ParmParse("xlo").addarr(
        "temperature", amrex::Vector<amrex::Real>{300.0});
    auto& temperature = sim().repo().declare_field("temperature", 1, 1, 2);
    kynema_sgf::BCScalar bc(temperature);
    bc();
    testing::internal::CaptureStdout();
    EXPECT_NO_THROW(validate());
    const auto output = testing::internal::GetCapturedStdout();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        const auto warning =
            output.find("Warning: xlo.temperature is overridden");
        ASSERT_NE(warning, std::string::npos);
        EXPECT_EQ(output.find("Warning:", warning + 1), std::string::npos);
    }
}

TEST_F(BCInputTest, modulated_power_law_supplies_velocity_and_temperature_only)
{
    initialize_mesh();
    auto& velocity = sim().repo().declare_field("velocity", 3, 1);
    kynema_sgf::BCVelocity bc(velocity);
    bc();
    scalar();
    amrex::ParmParse("ABL").addarr(
        "temperature_heights", amrex::Vector<amrex::Real>{0.0, 1.0});
    amrex::ParmParse("ABL").addarr(
        "temperature_values", amrex::Vector<amrex::Real>{300.0, 300.0});
    sim().field_boundary_manager().create("ModulatedPowerLaw", sim());
    EXPECT_NO_THROW(validate());
    scalar("density");
    EXPECT_THROW(validate(), amrex::RuntimeError);
}

} // namespace kynema_sgf_tests
