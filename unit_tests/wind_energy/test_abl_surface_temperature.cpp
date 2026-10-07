#include <fstream>
#include "abl_test_utils.H"
#include "src/wind_energy/ABLWallFunction.H"
#include "src/wind_energy/MOData.H"
#include "src/utilities/FieldPlaneAveraging.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {
void write_time_table(
    const std::string& fname,
    const std::string& header,
    const amrex::Vector<amrex::Real>& times,
    const amrex::Vector<amrex::Real>& values)
{
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(fname);
        os << header << '\n';
        for (amrex::Long n = 0; n < times.size(); ++n) {
            os << times[n] << '\t' << values[n] << '\n';
        }
    }
    amrex::ParallelDescriptor::Barrier();
}

//! Temperature of the Monin-Obukhov profile at height z
amrex::Real mo_temperature(const kynema_sgf::MOData& mo, const amrex::Real z)
{
    return mo.surf_temp -
           (mo.surf_temp_flux * mo.phi_h(z) / (mo.utau * mo.kappa));
}
} // namespace

class ABLSurfaceTemperatureTest : public ABLMeshTest
{
protected:
    void populate_parameters() override
    {
        ABLMeshTest::populate_parameters();
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", 2.5_rt);
    }

    //! Uniform 8 m/s wind at 300 K
    void setup_fields()
    {
        initialize_mesh();
        auto& pde_mgr = sim().pde_manager();
        pde_mgr.register_icns();
        pde_mgr.register_transport_pde("Temperature");
        auto& velocity = sim().repo().get_field("velocity");
        velocity.setVal(0.0_rt);
        velocity.setVal(8.0_rt, 0, 1);
        sim().repo().get_field("temperature").setVal(300.0_rt);
    }

    //! Move the simulation time to fixed_dt
    void advance_one_step()
    {
        auto& time = sim().time();
        time.new_timestep();
        time.set_current_cfl(0.0_rt, 0.0_rt, 0.0_rt);
        time.advance_time();
        time.new_timestep();
    }
};

TEST_F(ABLSurfaceTemperatureTest, heat_flux_timetable)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    write_time_table(
        "abl_heat_flux_table.txt", "time\tflux", {0.0_rt, 10.0_rt},
        {0.05_rt, 0.15_rt});
    populate_parameters();
    {
        amrex::ParmParse pp("ABL");
        pp.add("surface_temp_flux_timetable", "abl_heat_flux_table.txt");
    }
    setup_fields();

    kynema_sgf::VelPlaneAveraging vpa(sim(), 2, -1);
    kynema_sgf::FieldPlaneAveraging tpa(
        sim().repo().get_field("temperature"), sim().time(), 2, -1);
    vpa();
    tpa();

    kynema_sgf::ABLWallFunction wall_func(sim());
    wall_func.init_log_law_height(0);
    EXPECT_EQ(
        wall_func.mo().alg_type, kynema_sgf::MOData::ThetaCalcType::HEAT_FLUX);

    wall_func.update_umean(vpa, tpa);
    EXPECT_NEAR(wall_func.mo().surf_temp_flux, 0.05_rt, tol);

    advance_one_step();
    EXPECT_NEAR(sim().time().current_time(), 2.5_rt, tol);
    wall_func.update_umean(vpa, tpa);
    EXPECT_NEAR(wall_func.mo().surf_temp_flux, 0.075_rt, tol);
    // Unstable surface layer: the surface is warmer than the air
    EXPECT_GT(wall_func.mo().surf_temp, 300.0_rt);
}

TEST_F(ABLSurfaceTemperatureTest, near_surface_temperature_timetable)
{
    constexpr amrex::Real tol = 1.0e-3_rt;
    constexpr amrex::Real near_height = 2.0_rt;
    write_time_table(
        "abl_near_surface_table.txt", "time\ttemperature", {0.0_rt, 10.0_rt},
        {301.0_rt, 299.0_rt});
    populate_parameters();
    {
        amrex::ParmParse pp("ABL");
        pp.add("near_surface_temp_timetable", "abl_near_surface_table.txt");
        pp.add("near_surface_height", near_height);
    }
    setup_fields();

    kynema_sgf::VelPlaneAveraging vpa(sim(), 2, -1);
    kynema_sgf::FieldPlaneAveraging tpa(
        sim().repo().get_field("temperature"), sim().time(), 2, -1);
    vpa();
    tpa();

    kynema_sgf::ABLWallFunction wall_func(sim());
    wall_func.init_log_law_height(0);
    const auto& mo = wall_func.mo();
    EXPECT_EQ(
        mo.alg_type,
        kynema_sgf::MOData::ThetaCalcType::NEAR_SURFACE_TEMPERATURE);

    // Warm near-surface air under 300 K air: upward heat flux
    wall_func.update_umean(vpa, tpa);
    EXPECT_NEAR(mo.near_surf_temp, 301.0_rt, tol);
    EXPECT_GT(mo.surf_temp_flux, 0.0_rt);
    EXPECT_NEAR(mo_temperature(mo, near_height), 301.0_rt, tol);
    EXPECT_NEAR(mo_temperature(mo, mo.zref), 300.0_rt, tol);

    // Halfway through the table the air near the surface is at 300.5 K
    advance_one_step();
    wall_func.update_umean(vpa, tpa);
    EXPECT_NEAR(mo.near_surf_temp, 300.5_rt, tol);
    EXPECT_NEAR(mo_temperature(mo, near_height), 300.5_rt, tol);
    EXPECT_NEAR(mo_temperature(mo, mo.zref), 300.0_rt, tol);
}

TEST_F(ABLSurfaceTemperatureTest, near_surface_temperature_guards)
{
    write_time_table(
        "abl_near_surface_guard.txt", "time\ttemperature", {0.0_rt, 10.0_rt},
        {301.0_rt, 299.0_rt});
    populate_parameters();
    setup_fields();
    {
        // Same height as the log-law height
        amrex::ParmParse pp("ABL");
        pp.add("near_surface_temp_timetable", "abl_near_surface_guard.txt");
        pp.add("near_surface_height", 5.0_rt);
        pp.add("log_law_height", 5.0_rt);
        kynema_sgf::ABLWallFunction wall_func(sim());
        EXPECT_THROW(wall_func.init_log_law_height(0), std::runtime_error);
    }
    {
        // Combined with another surface temperature boundary condition
        amrex::ParmParse pp("ABL");
        pp.add("surface_temp_flux", 0.05_rt);
        EXPECT_THROW(
            kynema_sgf::ABLWallFunction wall_func(sim()), std::runtime_error);
    }
}

TEST_F(ABLSurfaceTemperatureTest, heat_flux_timetable_header_only)
{
    write_time_table("abl_header_only.txt", "time\tflux", {}, {});
    populate_parameters();
    {
        amrex::ParmParse pp("ABL");
        pp.add("surface_temp_flux_timetable", "abl_header_only.txt");
    }
    setup_fields();
    EXPECT_THROW(
        kynema_sgf::ABLWallFunction wall_func(sim()), std::runtime_error);
}

} // namespace kynema_sgf_tests
