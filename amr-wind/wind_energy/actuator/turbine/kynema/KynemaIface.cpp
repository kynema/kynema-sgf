#include <algorithm>
#include <cmath>
#include "amr-wind/wind_energy/actuator/turbine/external/ExtTurbIface.H"
#include "amr-wind/wind_energy/actuator/turbine/kynema/kynema_types.H"
#include "amr-wind/wind_energy/actuator/turbine/kynema/kynema_wrapper.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/SimTime.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/ncutils/nc_interface.H"

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_FileSystem.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace exw_kynema {
void build_turbine(
    kynema::interfaces::TurbineInterfaceBuilder& builder,
    const YAML::Node wio,
    const int n_blades,
    const int n_blade_nodes,
    const int n_tower_nodes,
    const amrex::Real rotor_speed_init,
    const amrex::Real generator_power_init,
    const amrex::Real wind_speed_init,
    const amrex::Real yaw_init,
    const amrex::Real generator_efficiency)
{
    // WindIO components
    const auto& wio_blade = wio["components"]["blade"];
    const auto& wio_tower = wio["components"]["tower"];
    const auto& wio_drivetrain = wio["components"]["drivetrain"];
    const auto& wio_hub = wio["components"]["hub"];
    const auto& wio_yaw = wio["components"]["yaw"];

    //--------------------------------------------------------------------------
    // Build Turbine
    //--------------------------------------------------------------------------

    // Get turbine builder
    auto& turbine_builder = builder.Turbine();
    turbine_builder.SetAzimuthAngle(0.)
        .SetRotorApexToHub(0.)
        .SetHubDiameter(wio_hub["diameter"].as<amrex::Real>())
        .SetConeAngle(
            wio_hub["cone_angle"].as<amrex::Real>() *
            std::numbers::pi_v<amrex::Real> / 180.0_rt)
        .SetShaftTiltAngle(
            wio_drivetrain["outer_shape"]["uptilt"].as<amrex::Real>() *
            std::numbers::pi_v<amrex::Real> / 180.0_rt)
        .SetTowerAxisToRotorApex(
            wio_drivetrain["outer_shape"]["overhang"].as<amrex::Real>())
        .SetTowerTopToRotorApex(
            wio_drivetrain["outer_shape"]["distance_tt_hub"].as<amrex::Real>())
        .SetGearBoxRatio(
            wio_drivetrain["gearbox"]["gear_ratio"].as<amrex::Real>())
        .SetRotorSpeed(rotor_speed_init)
        .SetGeneratorPower(generator_power_init)
        .SetHubWindSpeed(wind_speed_init)
        .SetNacelleYawAngle(yaw_init)
        .SetGeneratorEfficiency(generator_efficiency);

    //--------------------------------------------------------------------------
    // Build Blades
    //--------------------------------------------------------------------------

    // Loop through blades and set parameters
    for (auto j : std::views::iota(0U, static_cast<uint>(n_blades))) {
        // Get the blade builder
        auto& blade_builder = turbine_builder.Blade(j);

        // Set blade parameters
        blade_builder.SetElementOrder(n_blade_nodes - 1)
            .PrescribedRootMotion(false)
            .SetSectionRefinement(2);

        // Add reference axis coordinates (WindIO uses Z-axis as reference axis)
        const auto ref_axis = wio_blade["reference_axis"];
        const auto axis_grid =
            ref_axis["x"]["grid"].as<std::vector<amrex::Real>>();
        const auto x_values =
            ref_axis["x"]["values"].as<std::vector<amrex::Real>>();
        const auto y_values =
            ref_axis["y"]["values"].as<std::vector<amrex::Real>>();
        const auto z_values =
            ref_axis["z"]["values"].as<std::vector<amrex::Real>>();
        for (auto i : std::views::iota(0U, axis_grid.size())) {
            blade_builder.AddRefAxisPoint(
                axis_grid[i], {x_values[i], y_values[i], z_values[i]},
                kynema::interfaces::components::ReferenceAxisOrientation::Z);
        }

        // Add reference axis twist
        const auto blade_twist = wio_blade["outer_shape"]["twist"];
        const auto twist_grid =
            blade_twist["grid"].as<std::vector<amrex::Real>>();
        const auto twist_values =
            blade_twist["values"].as<std::vector<amrex::Real>>();
        for (auto i : std::views::iota(0U, twist_grid.size())) {
            blade_builder.AddRefAxisTwist(
                twist_grid[i],
                -twist_values[i] * std::numbers::pi_v<amrex::Real> / 180.0_rt);
        }

        const auto inertia_matrix =
            wio_blade["structure"]["elastic_properties"]["inertia_matrix"];
        const auto stiffness_matrix =
            wio_blade["structure"]["elastic_properties"]["stiffness_matrix"];

        // Add blade section properties
        const auto k_grid =
            stiffness_matrix["grid"].as<std::vector<amrex::Real>>();
        const auto m_grid =
            inertia_matrix["grid"].as<std::vector<amrex::Real>>();
        const auto n_sections = k_grid.size();
        if (m_grid.size() != k_grid.size()) {
            throw std::runtime_error(
                "stiffness and mass matrices not on same grid");
        }
        for (auto i : std::views::iota(0U, n_sections)) {
            if (std::abs(m_grid[i] - k_grid[i]) > 1.0e-8_rt) {
                throw std::runtime_error(
                    "stiffness and mass matrices not on same grid");
            }
            const auto mass = inertia_matrix["mass"][i].as<amrex::Real>();
            const auto cm_x = inertia_matrix["cm_x"][i].as<amrex::Real>();
            const auto cm_y = inertia_matrix["cm_y"][i].as<amrex::Real>();
            const auto i_cp = inertia_matrix["i_cp"][i].as<amrex::Real>();
            const auto i_edge = inertia_matrix["i_edge"][i].as<amrex::Real>();
            const auto i_flap = inertia_matrix["i_flap"][i].as<amrex::Real>();
            const auto i_plr = inertia_matrix["i_plr"][i].as<amrex::Real>();

            const auto k11 = stiffness_matrix["K11"][i].as<amrex::Real>();
            const auto k12 = stiffness_matrix["K12"][i].as<amrex::Real>();
            const auto k13 = stiffness_matrix["K13"][i].as<amrex::Real>();
            const auto k14 = stiffness_matrix["K14"][i].as<amrex::Real>();
            const auto k15 = stiffness_matrix["K15"][i].as<amrex::Real>();
            const auto k16 = stiffness_matrix["K16"][i].as<amrex::Real>();
            const auto k22 = stiffness_matrix["K22"][i].as<amrex::Real>();
            const auto k23 = stiffness_matrix["K23"][i].as<amrex::Real>();
            const auto k24 = stiffness_matrix["K24"][i].as<amrex::Real>();
            const auto k25 = stiffness_matrix["K25"][i].as<amrex::Real>();
            const auto k26 = stiffness_matrix["K26"][i].as<amrex::Real>();
            const auto k33 = stiffness_matrix["K33"][i].as<amrex::Real>();
            const auto k34 = stiffness_matrix["K34"][i].as<amrex::Real>();
            const auto k35 = stiffness_matrix["K35"][i].as<amrex::Real>();
            const auto k36 = stiffness_matrix["K36"][i].as<amrex::Real>();
            const auto k44 = stiffness_matrix["K44"][i].as<amrex::Real>();
            const auto k45 = stiffness_matrix["K45"][i].as<amrex::Real>();
            const auto k46 = stiffness_matrix["K46"][i].as<amrex::Real>();
            const auto k55 = stiffness_matrix["K55"][i].as<amrex::Real>();
            const auto k56 = stiffness_matrix["K56"][i].as<amrex::Real>();
            const auto k66 = stiffness_matrix["K66"][i].as<amrex::Real>();

            blade_builder.AddSection(
                m_grid[i],
                {{
                    {mass, 0., 0., 0., 0., -mass * cm_y},
                    {0., mass, 0., 0., 0., mass * cm_x},
                    {0., 0., mass, mass * cm_y, -mass * cm_x, 0.},
                    {0., 0., mass * cm_y, i_edge, -i_cp, 0.},
                    {0., 0., -mass * cm_x, -i_cp, i_flap, 0.},
                    {-mass * cm_y, mass * cm_x, 0., 0., 0., i_plr},
                }},
                {{
                    {k11, k12, k13, k14, k15, k16},
                    {k12, k22, k23, k24, k25, k26},
                    {k13, k23, k33, k34, k35, k36},
                    {k14, k24, k34, k44, k45, k46},
                    {k15, k25, k35, k45, k55, k56},
                    {k16, k26, k36, k46, k56, k66},
                }},
                kynema::interfaces::components::ReferenceAxisOrientation::Z);
        }
    }

    //--------------------------------------------------------------------------
    // Build Tower
    //--------------------------------------------------------------------------

    // Get the tower builder
    auto& tower_builder = turbine_builder.Tower();

    // Set tower parameters
    tower_builder
        .SetElementOrder(
            n_tower_nodes - 1) // Set element order to num nodes - 1
        .SetSectionRefinement(2)
        .PrescribedRootMotion(false); // Fix displacement of tower base node

    // Add reference axis coordinates (WindIO uses Z-axis as reference axis)
    const auto t_ref_axis = wio_tower["reference_axis"];
    const auto axis_grid =
        t_ref_axis["x"]["grid"].as<std::vector<amrex::Real>>();
    const auto x_values =
        t_ref_axis["x"]["values"].as<std::vector<amrex::Real>>();
    const auto y_values =
        t_ref_axis["y"]["values"].as<std::vector<amrex::Real>>();
    const auto z_values =
        t_ref_axis["z"]["values"].as<std::vector<amrex::Real>>();
    for (auto i : std::views::iota(0U, axis_grid.size())) {
        tower_builder.AddRefAxisPoint(
            axis_grid[i], {x_values[i], y_values[i], z_values[i]},
            kynema::interfaces::components::ReferenceAxisOrientation::Z);
    }

    // Set tower base position from first reference axis point
    turbine_builder.SetTowerBasePosition(
        {x_values[0], y_values[0], z_values[0], 1., 0., 0., 0.});

    // Add reference axis twist (zero for tower)
    tower_builder.AddRefAxisTwist(0.0_rt, 0.0_rt)
        .AddRefAxisTwist(1.0_rt, 0.0_rt);

    // Find the tower material properties
    const auto tower_wall_thickness =
        wio_tower["structure"]["layers"][0]["thickness"];
    const auto tower_material_name =
        wio_tower["structure"]["layers"][0]["material"].as<std::string>();

    YAML::Node tower_material;
    for (const auto& m : wio["materials"]) {
        if (m["name"] && m["name"].as<std::string>() == tower_material_name) {
            tower_material = m.as<YAML::Node>();
            break;
        }
    }
    if (!tower_material) {
        throw std::runtime_error(
            "Material '" + tower_material_name +
            "' not found in materials section");
    }

    // Add tower section properties
    const auto elastic_modulus = tower_material["E"].as<amrex::Real>();
    const auto shear_modulus = tower_material["G"].as<amrex::Real>();
    const auto poisson_ratio = tower_material["nu"].as<amrex::Real>();
    const auto density = tower_material["rho"].as<amrex::Real>();
    const auto tower_diameter = wio_tower["outer_shape"]["outer_diameter"];
    const auto tower_diameter_grid =
        tower_diameter["grid"].as<std::vector<amrex::Real>>();
    const auto tower_diameter_values =
        tower_diameter["values"].as<std::vector<amrex::Real>>();
    const auto tower_wall_thickness_values =
        tower_wall_thickness["values"].as<std::vector<amrex::Real>>();

    for (auto i : std::views::iota(0U, tower_diameter_grid.size())) {
        // Create section mass and stiffness matrices
        const auto section = kynema::beams::GenerateHollowCircleSection(
            tower_diameter_grid[i], elastic_modulus, shear_modulus, density,
            tower_diameter_values[i], tower_wall_thickness_values[i],
            poisson_ratio);

        // Add section
        tower_builder.AddSection(
            section.position, section.M_star, section.C_star,
            kynema::interfaces::components::ReferenceAxisOrientation::Z);
    }

    //--------------------------------------------------------------------------
    // Add mass elements
    //--------------------------------------------------------------------------

    // Get nacelle mass properties from WindIO
    const auto nacelle_props = wio_drivetrain["elastic_properties"];

    // Set the nacelle inertia matrix in the turbine builder
    if (nacelle_props) {
        // Get the nacelle mass, inertia properties, and location relative to
        // the tower top
        const auto nacelle_mass = nacelle_props["mass"].as<amrex::Real>();
        const auto nacelle_inertia =
            nacelle_props["inertia"].as<std::vector<amrex::Real>>();
        const auto nacelle_cm_offset =
            nacelle_props["location"].as<std::vector<amrex::Real>>();

        turbine_builder.SetNacelleInertiaMatrix(
            {{{nacelle_mass, 0., 0., 0., 0., 0.},
              {0., nacelle_mass, 0., 0., 0., 0.},
              {0., 0., nacelle_mass, 0., 0., 0.},
              {0., 0., 0., nacelle_inertia[0], nacelle_inertia[3],
               nacelle_inertia[4]},
              {0., 0., 0., nacelle_inertia[3], nacelle_inertia[1],
               nacelle_inertia[5]},
              {0., 0., 0., nacelle_inertia[4], nacelle_inertia[5],
               nacelle_inertia[2]}}},
            {nacelle_cm_offset[0], nacelle_cm_offset[1], nacelle_cm_offset[2]});
    }

    // Set yaw bearing inertia if yaw component exists
    if (wio_yaw) {
        // Get yaw bearing mass properties from WindIO
        const auto yaw_bearing_mass =
            wio_yaw["elastic_properties"]["mass"].as<amrex::Real>();
        const auto yaw_bearing_inertia =
            wio_yaw["elastic_properties"]["inertia"]
                .as<std::vector<amrex::Real>>();

        // Set the yaw bearing inertia matrix in the turbine builder
        turbine_builder.SetYawBearingInertiaMatrix(
            {{{yaw_bearing_mass, 0., 0., 0., 0., 0.},
              {0., yaw_bearing_mass, 0., 0., 0., 0.},
              {0., 0., yaw_bearing_mass, 0., 0., 0.},
              {0., 0., 0., yaw_bearing_inertia[0], 0., 0.},
              {0., 0., 0., 0., yaw_bearing_inertia[1], 0.},
              {0., 0., 0., 0., 0., yaw_bearing_inertia[2]}}});
    }

    // Get generator rotational inertia and gearbox ratio from WindIO
    const auto generator_inertia =
        wio_drivetrain["generator"]["elastic_properties"]["inertia"]
            .as<std::vector<amrex::Real>>();
    const auto gearbox_ratio =
        wio_drivetrain["gearbox"]["gear_ratio"].as<amrex::Real>();

    // Get hub mass properties from WindIO
    const auto hub_mass =
        wio_hub["elastic_properties"]["mass"].as<amrex::Real>();
    const auto hub_inertia =
        wio_hub["elastic_properties"]["inertia"].as<std::vector<amrex::Real>>();

    // Set the hub inertia matrix in the turbine builder
    turbine_builder.SetHubInertiaMatrix(
        {{{hub_mass, 0., 0., 0., 0., 0.},
          {0., hub_mass, 0., 0., 0., 0.},
          {0., 0., hub_mass, 0., 0., 0.},
          {0., 0., 0.,
           hub_inertia[0] +
               generator_inertia[0] * gearbox_ratio * gearbox_ratio,
           hub_inertia[3], hub_inertia[4]},
          {0., 0., 0., hub_inertia[3], hub_inertia[1], hub_inertia[5]},
          {0., 0., 0., hub_inertia[4], hub_inertia[5], hub_inertia[2]}}});
}

int build_aero(
    kynema::interfaces::TurbineInterfaceBuilder& builder, const YAML::Node wio)
{
    //--------------------------------------------------------------------------
    // Build Aerodynamics
    //--------------------------------------------------------------------------

    auto& aero_builder = builder.Aerodynamics()
                             .EnableAero()
                             .SetNumberOfAirfoils(1UL)
                             .SetAirfoilToBladeMap(std::array{0UL, 0UL, 0UL});

    const auto& airfoil_io = wio["airfoils"];
    auto aero_sections =
        std::vector<kynema::interfaces::components::AerodynamicSection>{};
    auto id = 0UL;
    for (const auto& af : airfoil_io) {
        const auto s = af["spanwise_position"].as<amrex::Real>();
        const auto chord = af["chord"].as<amrex::Real>();
        const auto twist = af["twist"].as<amrex::Real>() *
                           std::numbers::pi_v<amrex::Real> / 180.0_rt;
        const auto section_offset_x = af["section_offset_x"].as<amrex::Real>();
        const auto section_offset_y = af["section_offset_y"].as<amrex::Real>();
        const auto aerodynamic_center =
            af["aerodynamic_center"].as<amrex::Real>();
        auto aoa = af["polars"][0]["re_sets"][0]["cl"]["grid"]
                       .as<std::vector<amrex::Real>>();
        std::ranges::transform(aoa, std::begin(aoa), [](auto degrees) {
            return degrees * std::numbers::pi_v<amrex::Real> / 180.0_rt;
        });
        const auto cl = af["polars"][0]["re_sets"][0]["cl"]["values"]
                            .as<std::vector<amrex::Real>>();
        const auto cd = af["polars"][0]["re_sets"][0]["cd"]["values"]
                            .as<std::vector<amrex::Real>>();
        const auto cm = af["polars"][0]["re_sets"][0]["cm"]["values"]
                            .as<std::vector<amrex::Real>>();
        aero_sections.emplace_back(
            id, s, chord, section_offset_x, section_offset_y,
            aerodynamic_center, twist, aoa, cl, cd, cm);
        ++id;
    }

    aero_builder.SetAirfoilSections(0UL, aero_sections);

    return (int)aero_sections.size();
}

void build_controller(
    kynema::interfaces::TurbineInterfaceBuilder& builder,
    const std::string controller_shared_lib_path,
    const std::string controller_input_file,
    const std::string turbine_label)
{
    const auto controller_function_name = std::string{"DISCON"};
    const auto controller_output_file = "controller_" + turbine_label;

    auto controller_builder = builder.Controller()
                                  .SetLibraryPath(controller_shared_lib_path)
                                  .SetFunctionName(controller_function_name)
                                  .SetInputFilePath(controller_input_file)
                                  .SetControllerInput(controller_output_file);
}

void update_turbine(::ext_turb::KynemaTurbine& fi, bool advance)
{
    ++fi.substep_counter;
    // Operations that only need to be done once
    if (fi.substep_counter == 1) {
        fi.interface->Aerodynamics().CalculateMotion(
            fi.interface->GetHostState());
        // copy fluid velocity to turbine solver (set inflow)
        fi.pass_fluid_velocity_and_hub_load();
        fi.interface->Aerodynamics().CalculateAerodynamicLoads(
            fi.fluid_density);
        fi.interface->Aerodynamics().CalculateNodalLoads();
    }
    if (advance) {
        // individual turbine step, do not output every step
        bool converged = fi.interface->Step();
        if (!converged) {
            amrex::Abort("Kynema did not converge\n");
        }
        if (fi.controller_input_file.size() > 0) {
            const amrex::Real t = fi.time_index * fi.dt_ext;
            fi.interface->ApplyController(t);
        }
    }

    if (fi.substep_counter == fi.num_substeps) {
        fi.interface->Aerodynamics().CalculateMotion(
            fi.interface->GetHostState());
        fi.substep_counter = 0;
    } else if (!advance) {
        fi.substep_counter = 0;
    }

    if (fi.substep_counter == 0) {
        // Output once per amr-wind timestep
        fi.interface->OpenOutputFile();
        fi.interface->WriteOutput();
        fi.interface->CloseOutputFile();
        // Populate buffers with turbine data
        fi.populate_buffers();
    }
}
} // namespace exw_kynema

namespace ext_turb {

template <>
ExtTurbIface<KynemaTurbine, KynemaSolverData>::~ExtTurbIface()
{
    //! Do deallocation if necessary
}

// !! This doesn't get used by the actual code, just by unit tests !! //
template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::parse_inputs(
    const amr_wind::CFDSim& sim, const std::string& inp_name)
{
    amrex::ParmParse pp(inp_name);

    const auto& time = sim.time();
    //! Check that the user has not enabled adaptive timestepping
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        !time.adaptive_timestep(),
        "Adaptive time-stepping not supported when Kynema is enabled");

    m_dt_cfd = time.delta_t();
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::allocate_ext_turbines()
{
    BL_PROFILE("amr-wind::KynemaIface::allocate_turbines");

    m_is_initialized = true;

    // Get solver parameters
    {
        amrex::ParmParse pp("Kynema");
        pp.query("damping_factor", m_solver_data.damping_factor);
        pp.query("max_nonlinear_iterations", m_solver_data.nl_iter_max);
        pp.query("abs_err_tol", m_solver_data.abs_err_tol);
        pp.query("rel_err_tol", m_solver_data.rel_err_tol);
    }

    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_solver_data.gravity);
    }
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::init_solution(
    const int local_id)
{
    BL_PROFILE("amr-wind::KynemaIface::init_solution");
    AMREX_ALWAYS_ASSERT(local_id < static_cast<int>(m_turbine_data.size()));
    AMREX_ALWAYS_ASSERT(m_is_initialized);

    auto& fi = *m_turbine_data[local_id];

    ::exw_kynema::update_turbine(fi, false);

    fi.is_solution0 = false;
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::get_hub_stats(
    const int local_id)
{
    BL_PROFILE("amr-wind::KynemaIface::get_hub_stats");

    auto& fi = *m_turbine_data[local_id];

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        fi.hub_abs_pos[dir] = fi.interface->Turbine().hub_node.position[dir];
        fi.hub_rot_vel[dir] =
            fi.interface->Turbine().hub_node.velocity[dir + AMREX_SPACEDIM];
    }
    // Orientation of hub is already in turbine data as first point
    for (int comp = 0; comp < AMREX_SPACEDIM * AMREX_SPACEDIM; ++comp) {
        fi.hub_orient[comp] = fi.orientation()[comp];
    }
}

#ifdef AMR_WIND_USE_KYNEMA

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::prepare_netcdf_file(
    KynemaTurbine& fi)
{
#ifdef AMR_WIND_USE_NETCDF
    BL_PROFILE("amr-wind::KynemaIface::prepare_netcdf_file");
    if (!amrex::UtilCreateDirectory(m_solver_data.output_dir, 0755)) {
        amrex::CreateDirectoryFailed(m_solver_data.output_dir);
    }

    const std::string fname =
        m_solver_data.output_dir + "/" + fi.tlabel + ".nc";

    // Don't overwrite existing
    if (amrex::FileSystem::Exists(fname)) {
        return;
    }

    auto ncf = ncutils::NCFile::create(fname, NC_CLOBBER | NC_NETCDF4);
    const std::string nt_name = "num_time_steps";
    const std::string np_name = "num_vel_points";
    ncf.enter_def_mode();
    ncf.put_attr("title", "AMR-Wind Kynema velocity data");
    ncf.put_attr("AMR-Wind_version", amr_wind::ioutils::amr_wind_version());
    ncf.put_attr("created_on", amr_wind::ioutils::timestamp());
    ncf.def_dim(nt_name, NC_UNLIMITED);
    ncf.def_dim(np_name, fi.length_fluid_velocity(0));
    ncf.def_dim("ndim", AMREX_SPACEDIM);
    ncf.def_var("time", NC_FLOAT, {nt_name});
    ncf.def_var("xco", NC_FLOAT, {np_name});
    ncf.def_var("yco", NC_FLOAT, {np_name});
    ncf.def_var("zco", NC_FLOAT, {np_name});
    ncf.def_var("uvel", NC_FLOAT, {nt_name, np_name});
    ncf.def_var("vvel", NC_FLOAT, {nt_name, np_name});
    ncf.def_var("wvel", NC_FLOAT, {nt_name, np_name});
    ncf.exit_def_mode();

    {
        const auto npts = static_cast<size_t>(fi.length_fluid_velocity(0));
        auto xco = ncf.var("xco");
        xco.put(fi.position_at_vel(0), {0}, {npts});
        auto yco = ncf.var("yco");
        yco.put(fi.position_at_vel(1), {0}, {npts});
        auto zco = ncf.var("zco");
        zco.put(fi.position_at_vel(2), {0}, {npts});
    }
#else
    amrex::ignore_unused(fi);
    amrex::OutStream() << "WARNING: KynemaIface: NetCDF support was not "
                          "enabled during compile "
                          "time. KynemaIface cannot support restart."
                       << std::endl;
#endif
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::write_velocity_data(
    const KynemaTurbine& fi)
{
#ifdef AMR_WIND_USE_NETCDF
    BL_PROFILE("amr-wind::KynemaIface::write_velocity_data");
    const std::string fname =
        m_solver_data.output_dir + "/" + fi.tlabel + ".nc";
    auto ncf = ncutils::NCFile::open(fname, NC_WRITE);
    const std::string nt_name = "num_time_steps";
    const size_t nt = ncf.dim(nt_name).len();
    const auto npts = static_cast<size_t>(fi.length_fluid_velocity(0));

    const amrex::Real time = fi.time_index * fi.dt_ext;
    ncf.var("time").put(&time, {nt}, {1});
    ncf.var("uvel").put(fi.fluid_velocity(0), {nt, 0}, {1, npts});
    ncf.var("vvel").put(fi.fluid_velocity(1), {nt, 0}, {1, npts});
    ncf.var("wvel").put(fi.fluid_velocity(2), {nt, 0}, {1, npts});
#else
    amrex::ignore_unused(fi);
#endif
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::read_velocity_data(
    KynemaTurbine& fi, const ncutils::NCFile& ncf, const size_t tid)
{
#ifdef AMR_WIND_USE_NETCDF
    const auto nt = static_cast<size_t>(tid);
    const auto npts = static_cast<size_t>(fi.length_fluid_velocity(0));

    ncf.var("uvel").get(fi.fluid_velocity(0), {nt, 0}, {1, npts});
    ncf.var("vvel").get(fi.fluid_velocity(1), {nt, 0}, {1, npts});
    ncf.var("wvel").get(fi.fluid_velocity(2), {nt, 0}, {1, npts});
#else
    amrex::ignore_unused(fi);
    amrex::Abort(
        "KynemaIface::read_velocity_data: AMR-Wind was not compiled with "
        "NetCDF "
        "support");
#endif
}

#else

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::prepare_netcdf_file(
    KynemaTurbine& /*unused*/)
{}
template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::write_velocity_data(
    const KynemaTurbine& /*unused*/)
{}
template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::read_velocity_data(
    KynemaTurbine& /*unused*/,
    const ncutils::NCFile& /*unused*/,
    const size_t /*unused*/)
{}

#endif

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::do_turbine_step(
    KynemaTurbine& fi)
{
    ::exw_kynema::update_turbine(fi, true);
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::write_turbine_checkpoint(
    int& tid)
{
    // write checkpoint
}

template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::ext_init_turbine(
    KynemaTurbine& fi)
{
    BL_PROFILE("amr-wind::KynemaIface::init_turbine");

    auto builder = kynema::interfaces::TurbineInterfaceBuilder{};

    builder.Solution()
        .EnableDynamicSolve()
        .SetTimeStep(fi.dt_ext)
        .SetDampingFactor(m_solver_data.damping_factor)
        .SetGravity(
            {m_solver_data.gravity[0], m_solver_data.gravity[1],
             m_solver_data.gravity[2]})
        .SetMaximumNonlinearIterations(
            static_cast<size_t>(m_solver_data.nl_iter_max))
        .SetAbsoluteErrorTolerance(m_solver_data.abs_err_tol)
        .SetRelativeErrorTolerance(m_solver_data.rel_err_tol);

    const YAML::Node wio = YAML::LoadFile(fi.input_file);

    constexpr int num_pts_tower_struct{11};

    // Builds turbine, including blades, nacelle, and tower
    exw_kynema::build_turbine(
        builder, wio, fi.num_blades, fi.num_blade_elem, num_pts_tower_struct,
        fi.rotational_speed, fi.generator_power, fi.wind_speed, fi.yaw,
        fi.generator_efficiency);

    auto n_aero_sections = exw_kynema::build_aero(builder, wio);

    if (n_aero_sections != fi.num_pts_blade) {
        amrex::Abort(
            "KynemaIface: number of points per blade (from AMR-Wind input "
            "file) does not match number of aerodynamic sections per blade "
            "(from Kynema input file).");
    }

    if (fi.num_pts_tower != 0) {
        amrex::Abort(
            "KynemaIface: aerodynamic tower points not yet supported; number "
            "of tower points in input file must be 0.");
    }

    if (fi.controller_input_file.size() > 0) {
        exw_kynema::build_controller(
            builder, fi.controller_shared_lib_path, fi.controller_input_file,
            fi.tlabel);
    }

    // Create output
    builder.Outputs().SetOutputFilePath("kynema_" + fi.tlabel);

    fi.interface = std::make_unique<kynema::interfaces::TurbineInterface>(
        builder.Solution().Input(), builder.Turbine().Input(),
        builder.Aerodynamics().Input(), builder.Controller().Input(),
        builder.Outputs().Config());

    // Close file
    fi.interface->CloseOutputFile();

    // Determine the number of substeps for Kynema per CFD timestep
    fi.num_substeps = static_cast<int>(std::floor(fi.dt_cfd / fi.dt_ext));

    AMREX_ALWAYS_ASSERT(fi.num_substeps > 0);
    // Check that the time step sizes are consistent and Kynema advances at an
    // integral multiple of CFD timestep
    amrex::Real dt_err =
        fi.dt_cfd / (static_cast<amrex::Real>(fi.num_substeps) * fi.dt_ext) -
        1.0_rt;
    if (dt_err > 1.0e-12_rt) {
        amrex::Abort(
            "KynemaIFace: Kynema timestep is not an integral "
            "multiple of CFD timestep");
    }

    fi.allocate_buffers();
    ::exw_kynema::update_turbine(fi, false);
}

// cppcheck-suppress constParameterReference
// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::ext_replay_turbine(
    KynemaTurbine& fi)
{

    // Do we even do this???
}

// NOLINTNEXTLINE(readability-convert-member-functions-to-static)
template <>
void ExtTurbIface<KynemaTurbine, KynemaSolverData>::ext_restart_turbine(
    KynemaTurbine& fi)
{
    BL_PROFILE("amr-wind::KynemaIface::restart_turbine");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        amrex::FileSystem::Exists(fi.checkpoint_file + ".chkp"),
        "KynemaIface: Cannot find Kynema checkpoint file: " +
            fi.checkpoint_file);

    // Determine the number of substeps for Kynema per CFD timestep
    fi.num_substeps = static_cast<int>(std::floor(fi.dt_cfd / fi.dt_ext));

    AMREX_ALWAYS_ASSERT(fi.num_substeps > 0);
    // Check that the time step sizes are consistent and Kynema advances at an
    // integral multiple of CFD timestep
    amrex::Real dt_err =
        fi.dt_cfd / (static_cast<amrex::Real>(fi.num_substeps) * fi.dt_ext) -
        1.0_rt;
    if (dt_err > 1.0e-4_rt) {
        amrex::Abort(
            "KynemaIFace: Kynema timestep is not an integral "
            "multiple of CFD timestep");
    }
}

template class ExtTurbIface<KynemaTurbine, KynemaSolverData>;

} // namespace ext_turb
