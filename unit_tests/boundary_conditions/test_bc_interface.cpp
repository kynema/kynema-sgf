#include "ks_test_utils/MeshTest.H"
#include "src/boundary_conditions/BCInterface.H"
#include "src/boundary_conditions/field_boundary_fill/BoundaryPlane.H"
#include "src/boundary_conditions/field_boundary_fill/ModulatedPowerLaw.H"
#include "src/boundary_conditions/field_boundary_fill/OceanWavesBoundary.H"
#include "src/equation_systems/vof/vof_bciface.H"
#include "src/incflo.H"
#include "AMReX_Utility.H"

#include <filesystem>

namespace kynema_sgf_tests {

namespace {
const amrex::Orientation xlo(0, amrex::Orientation::low);
const amrex::Orientation xhi(0, amrex::Orientation::high);

template <typename Action>
void expect_missing_inflow(Action action, const std::string& key)
{
    try {
        action();
        FAIL() << "Expected missing inflow input " << key;
    } catch (const amrex::RuntimeError& error) {
        EXPECT_NE(
            std::string(error.what())
                .find("BC: missing required inflow value " + key),
            std::string::npos);
    }
}

amrex::Real ghost_value(
    const kynema_sgf::Field& field,
    const amrex::IntVect& iv,
    const int component = 0)
{
    auto value = amrex::ReduceSum(
        field(0), 1,
        [=] AMREX_GPU_HOST_DEVICE(
            const amrex::Box& bx,
            const amrex::Array4<const amrex::Real>& arr) -> amrex::Real {
            return bx.contains(iv) ? arr(iv, component) : 0.0;
        });
    amrex::ParallelDescriptor::ReduceRealSum(value);
    return value;
}
} // namespace

class BCInputTest : public MeshTest
{
public:
    void TearDown() override
    {
        m_mesh.reset();
        if (!m_plane_directory.empty() &&
            amrex::ParallelDescriptor::IOProcessor()) {
            std::error_code error;
            std::filesystem::remove_all(m_plane_directory, error);
            EXPECT_FALSE(error) << error.message();
        }
        MeshTest::TearDown();
    }

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

    void validate(
        const kynema_sgf::InflowValidationStage stage =
            kynema_sgf::InflowValidationStage::initialized)
    {
        kynema_sgf::BCIface::validate_inflow_values(
            sim().repo(), sim().field_boundary_manager(), stage);
    }

    void configure_static_plane(
        const amrex::Vector<std::string>& fields,
        const std::string& fluid_phase = "both")
    {
        std::string directory = "bc_input_test_" + amrex::UniqueString();
        int length = static_cast<int>(directory.size());
        const int root = amrex::ParallelDescriptor::IOProcessorNumber();
        amrex::ParallelDescriptor::Bcast(&length, 1, root);
        directory.resize(length);
        amrex::ParallelDescriptor::Bcast(directory.data(), length, root);
        if (amrex::ParallelDescriptor::IOProcessor()) {
            ASSERT_TRUE(std::filesystem::create_directory(directory));
        }
        m_plane_directory = directory;
        amrex::ParmParse pp("BoundaryPlane");
        pp.add("io_mode", 1);
        pp.add("file", directory + "/planes");
        pp.addarr("var_names", fields);
        pp.add("fluid_phase", fluid_phase);
        pp.add("output_format", std::string("native"));
        pp.add("output_and_use_initial_plane", true);
        pp.addarr("planes", amrex::Vector<std::string>{"xlo"});
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

    std::string m_plane_directory;
};

TEST_F(BCInputTest, driver_rejects_missing_velocity_before_meshing)
{
    populate_parameters();
    reset_prob_domain();
    amrex::ParmParse("xlo").addarr("density", amrex::Vector<amrex::Real>{1.0});
    expect_missing_inflow([]() { incflo driver; }, "xlo.velocity");
}

TEST_F(BCInputTest, driver_rejects_missing_density_before_meshing)
{
    populate_parameters();
    reset_prob_domain();
    amrex::ParmParse("xlo").addarr(
        "velocity", amrex::Vector<amrex::Real>{0.0, 0.0, 0.0});
    expect_missing_inflow([]() { incflo driver; }, "xlo.density");
}

TEST_F(BCInputTest, driver_rejects_missing_vof_before_meshing)
{
    populate_parameters();
    reset_prob_domain();
    amrex::ParmParse("incflo").add("use_godunov", 1);
    amrex::ParmParse("incflo").addarr(
        "physics", amrex::Vector<std::string>{"MultiPhase"});
    amrex::ParmParse("xlo").addarr(
        "velocity", amrex::Vector<amrex::Real>{0.0, 0.0, 0.0});
    amrex::ParmParse("xlo").addarr("density", amrex::Vector<amrex::Real>{1.0});
    expect_missing_inflow([]() { incflo driver; }, "xlo.vof");
}

TEST_F(BCInputTest, driver_checks_inflow_values_after_boundary_setup)
{
    populate_parameters();
    reset_prob_domain();
    amrex::ParmParse("xlo").addarr(
        "velocity", amrex::Vector<amrex::Real>{0.0, 0.0, 0.0});
    amrex::ParmParse("xlo").addarr("density", amrex::Vector<amrex::Real>{1.0});
    incflo driver;
    driver.init_mesh();
    auto& late_scalar = driver.repo().declare_field("late_scalar", 1, 1);
    kynema_sgf::BCScalar bc(late_scalar);
    bc();
    expect_missing_inflow(
        [&]() { driver.init_kynema_sgf_modules(); }, "xlo.late_scalar");
}

TEST_F(
    BCInputTest,
    driver_accepts_initial_plane_without_velocity_or_temperature_constants)
{
    populate_parameters();
    reset_prob_domain();
    amrex::ParmParse("incflo").addarr(
        "field_boundaries", amrex::Vector<std::string>{"BoundaryPlane"});
    amrex::ParmParse("xlo").addarr("density", amrex::Vector<amrex::Real>{1.0});
    configure_static_plane({"velocity", "temperature"});
    incflo driver;
    driver.sim().pde_manager().register_transport_pde("Temperature");
    driver.init_mesh();
    driver.velocity().setVal(2.0);
    driver.density().setVal(1.0);
    driver.temperature().setVal(300.0);
    EXPECT_NO_THROW(driver.init_kynema_sgf_modules());
    EXPECT_EQ(ghost_value(driver.velocity(), {-1, 0, 0}), 2.0);
    EXPECT_EQ(ghost_value(driver.temperature(), {-1, 0, 0}), 300.0);
}

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
    EXPECT_NO_THROW(validate(kynema_sgf::InflowValidationStage::configuration));
    EXPECT_THROW(validate(), amrex::RuntimeError);
}

TEST_F(BCInputTest, configuration_validation_defers_plane_data_loading)
{
    populate_parameters();
    create_mesh_instance();
    auto& temperature = scalar();
    amrex::ParmParse pp("BoundaryPlane");
    pp.add("io_mode", 1);
    pp.add("file", std::string("not_loaded"));
    pp.addarr("var_names", amrex::Vector<std::string>{"temperature"});
    auto& boundary =
        sim().field_boundary_manager().create("BoundaryPlane", sim());
    EXPECT_TRUE(boundary.may_provide_inflow_values(temperature, xlo));
    EXPECT_FALSE(boundary.provides_inflow_values(temperature, xlo));
    EXPECT_NO_THROW(validate(kynema_sgf::InflowValidationStage::configuration));
    scalar("density");
    EXPECT_THROW(
        validate(kynema_sgf::InflowValidationStage::configuration),
        amrex::RuntimeError);
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

TEST_F(BCInputTest, vof_and_multiphase_density_do_not_use_scalar_profiles)
{
    initialize_mesh();
    amrex::ParmParse("xlo").add(
        "vof.inflow_type", std::string("LinearProfile"));
    amrex::ParmParse("xlo").add(
        "density.inflow_type", std::string("LinearProfile"));
    auto& vof = sim().repo().declare_field("vof", 1, 1);
    kynema_sgf::BCVOF bc_vof(vof);
    bc_vof();
    expect_missing_inflow([&]() { validate(); }, "xlo.vof");
    amrex::ParmParse("xlo").addarr("vof", amrex::Vector<amrex::Real>{0.0});
    auto& density = sim().repo().declare_field("density", 1, 1);
    kynema_sgf::BCDensityMod bc_density(density);
    bc_density();
    expect_missing_inflow([&]() { validate(); }, "xlo.density");
}

TEST_F(BCInputTest, redundant_constants_produce_a_warning)
{
    initialize_mesh();
    amrex::ParmParse("xlo").addarr(
        "temperature", amrex::Vector<amrex::Real>{300.0});
    scalar();
    plane({"temperature"});
    testing::internal::CaptureStdout();
    EXPECT_NO_THROW(validate(kynema_sgf::InflowValidationStage::configuration));
    EXPECT_TRUE(testing::internal::GetCapturedStdout().empty());
    testing::internal::CaptureStdout();
    EXPECT_NO_THROW(validate());
    const auto output = testing::internal::GetCapturedStdout();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        EXPECT_NE(
            output.find("Warning: xlo.temperature is overridden"),
            std::string::npos);
        EXPECT_NE(output.find("not needed"), std::string::npos);
    }
}

TEST_F(BCInputTest, gas_plane_and_ocean_waves_masks_can_leave_a_gap)
{
    populate_parameters();
    amrex::ParmParse("incflo").add("use_godunov", 1);
    amrex::ParmParse("incflo").addarr(
        "physics", amrex::Vector<std::string>{"MultiPhase", "OceanWaves"});
    initialize_mesh();
    sim().pde_manager().register_icns();
    sim().init_physics();
    sim().init_field_boundaries();
    auto& repo = sim().repo();
    auto& velocity = repo.get_field("velocity");
    auto& waves =
        sim().field_boundary_manager().get<kynema_sgf::OceanWavesBoundary>();
    auto& vof = repo.get_field("vof");
    vof.setVal(1.0);
    repo.get_field("ow_vof").setVal(0.0);
    repo.get_field("ow_velocity").setVal(11.0);
    velocity.setVal(42.0);
    configure_static_plane({"velocity"}, "gas");
    auto& boundary =
        sim().field_boundary_manager().create("BoundaryPlane", sim());
    boundary.post_init_actions();

    velocity.setVal(7.0);
    boundary.set_velocity(0, 0.0, velocity, velocity(0));
    waves.set_velocity(0, 0.0, velocity, velocity(0));
    EXPECT_EQ(ghost_value(velocity, {-1, 0, 0}), 7.0);
    EXPECT_FALSE(boundary.provides_inflow_values(velocity, xlo));
    EXPECT_FALSE(waves.provides_inflow_values(velocity, xlo));
    EXPECT_THROW(validate(), amrex::RuntimeError);
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

    plane({"temperature"});
    testing::internal::CaptureStdout();
    EXPECT_NO_THROW(validate());
    const auto plane_output = testing::internal::GetCapturedStdout();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        EXPECT_NE(plane_output.find("an inflow profile"), std::string::npos);
        EXPECT_EQ(plane_output.find("fallback"), std::string::npos);
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
