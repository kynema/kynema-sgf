#include <array>
#include <cmath>
#include <numbers>
#include <ranges>
#include <stdexcept>
#include <string>
#include <vector>

// This file is compiled with the host compiler when Kokkos is host-only (see
// CMakeLists.txt), so it must not include AMReX headers that contain GPU code
// (e.g. AMReX_Vector.H, AMReX_Array.H); standard containers are used instead.
#include "src/wind_energy/actuator/turbine/kynema_fmb/kynema_fmb_host.H"
#include "src/wind_energy/actuator/turbine/kynema_fmb/kynema_fmb_wrapper.H"

using namespace amrex::literals;

namespace sgf_fmb {

void TurbineInterfaceDeleter::operator()(TurbineInterface* ptr) const
{
    delete ptr;
}

namespace {
//! Index of the aerodynamic body for actuator point id (hub is point 0)
int body_index(TurbineInterface& ti, int id)
{
    // Considering hub, which is 1 node, subtract
    id -= 1;
    int b_id = 0;
    int id_diff = id;
    const auto& bodies = ti.Aerodynamics().bodies;
    for (int b = 0; b < static_cast<int>(bodies.size()); ++b) {
        id_diff -= static_cast<int>(bodies[b].loads.size());
        if (id_diff < 0) {
            b_id = b;
            break;
        }
    }
    return b_id;
}

//! Index of the node within its aerodynamic body for actuator point id
int node_index(TurbineInterface& ti, int id)
{
    // Considering hub, which is 1 node, subtract
    id -= 1;
    int n_id = 0;
    int id_diff = id;
    const auto& bodies = ti.Aerodynamics().bodies;
    for (int b = 0; b < static_cast<int>(bodies.size()); ++b) {
        id_diff -= static_cast<int>(bodies[b].loads.size());
        if (id_diff < 0) {
            n_id = id_diff + static_cast<int>(bodies[b].loads.size());
            break;
        }
    }
    return n_id;
}
} // namespace

int num_force_points(TurbineInterface& ti)
{
    unsigned int n_sect{1};
    for (const auto& body : ti.Aerodynamics().bodies) {
        n_sect += body.loads.size();
    }
    return static_cast<int>(n_sect);
}

void populate_buffers(
    TurbineInterface& ti,
    const int npts,
    const bool one_way_coupled,
    float* chord,
    float* pos_force,
    float* pos_vel,
    float* solid_vel,
    float* force,
    float* orientation)
{
    const auto& hub_node = ti.Turbine().hub_node;
    const auto& bodies = ti.Aerodynamics().bodies;

    // Get hub information for nacelle force
    // Chord is not used, set it to zero
    chord[0] = 0.;
    // Force is calculated by Kynema-SGF, passed to Kynema-FMB
    // Get position directly from node information
    pos_force[0] = hub_node.position[0];
    pos_force[npts] = hub_node.position[1];
    pos_force[2 * npts] = hub_node.position[2];
    pos_vel[0] = pos_force[0];
    pos_vel[npts] = pos_force[npts];
    pos_vel[2 * npts] = pos_force[2 * npts];
    // Get solid velocity as well, not used at the moment though
    solid_vel[0] = hub_node.velocity[0];
    solid_vel[npts] = hub_node.velocity[1];
    solid_vel[2 * npts] = hub_node.velocity[2];

    // Get blade and tower information
    int idx = 0;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        // Skip hub node
        ++idx;
        for (int n = 1; n < npts; ++n) {
            const auto& body = bodies[body_index(ti, n)];
            const int in = node_index(ti, n);
            if (dir == 0) {
                chord[n] = (float)body.chord[in];
            }
            if (!one_way_coupled) {
                force[idx] = (float)body.loads[in][dir];
            }
            pos_force[idx] = (float)body.x_motion[in][dir];
            pos_vel[idx] = pos_force[idx];
            solid_vel[idx] = (float)body.v_motion[in][dir];
            ++idx;
        }
    }
    // Orientation stores a 3x3 tensor (9 elements) for every point
    std::array<float, 4> q;
    for (int n = 0; n < npts; ++n) {
        // convert quaternions to rotation matrix
        for (int iq = 0; iq < 4; ++iq) {
            if (n == 0) {
                q[iq] = (float)hub_node.position[iq + 3];
            } else {
                q[iq] = (float)bodies[body_index(ti, n)]
                            .qqr_motion_map[node_index(ti, n)][iq];
            }
        }
        idx = n * 9;
        orientation[idx + 0] =
            q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
        orientation[idx + 1] = 2. * (q[1] * q[2] - q[0] * q[3]);
        orientation[idx + 2] = 2. * (q[1] * q[3] + q[0] * q[2]);
        orientation[idx + 3] = 2. * (q[1] * q[2] + q[0] * q[3]);
        orientation[idx + 4] =
            (q[0] * q[0] - q[1] * q[1] + q[2] * q[2] - q[3] * q[3]);
        orientation[idx + 5] = 2. * (q[2] * q[3] - q[0] * q[1]);
        orientation[idx + 6] = 2. * (q[1] * q[3] - q[0] * q[2]);
        orientation[idx + 7] = 2. * (q[2] * q[3] + q[0] * q[1]);
        orientation[idx + 8] =
            (q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3]);
    }
}

void pass_fluid_velocity_and_hub_load(
    TurbineInterface& ti,
    const int npts,
    const float* fluid_vel,
    const float* force)
{
    // Fluid velocity at hub does not need to be passed, loads computed
    // internally in Kynema-SGF
    int idx = 0;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        // Skip hub node
        ++idx;
        for (int n = 1; n < npts; ++n) {
            ti.Aerodynamics()
                .bodies[body_index(ti, n)]
                .v_inflow[node_index(ti, n)][dir] = fluid_vel[idx];
            ++idx;
        }
    }
    // Hub load (also called nacelle_force) does need to be passed
    ti.Turbine().hub_node.loads[0] = force[0];
    ti.Turbine().hub_node.loads[1] = force[npts];
    ti.Turbine().hub_node.loads[2] = force[2 * npts];
    // Fluid velocity at hub is passed separately for the sake of controller
    std::array<amrex::Real, 3> hub_inflow{
        fluid_vel[0], fluid_vel[npts], fluid_vel[2 * npts]};
    ti.SetHubInflow(hub_inflow);
}

void save_and_zero_aerodynamic_loads(
    TurbineInterface& ti, const int npts, float* force)
{
    // Zero blade and tower aero loads
    int idx = 0;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        // Skip hub node
        ++idx;
        for (int n = 1; n < npts; ++n) {
            auto& load =
                ti.Aerodynamics().bodies[body_index(ti, n)].loads[node_index(
                    ti, n)][dir];
            force[idx] = (float)load;
            load = 0.0_rt;
            ++idx;
        }
    }
    // Hub load
    ti.Turbine().hub_node.loads[0] = 0.0_rt;
    ti.Turbine().hub_node.loads[1] = 0.0_rt;
    ti.Turbine().hub_node.loads[2] = 0.0_rt;
}

namespace {
void build_turbine(
    kynema_fmb::interfaces::TurbineInterfaceBuilder& builder,
    const YAML::Node wio,
    const int n_blades,
    const int n_blade_nodes,
    const int n_tower_nodes,
    const int section_refinement_blade,
    const int section_refinement_tower,
    const amrex::Real rotor_speed_init,
    const amrex::Real generator_power_init,
    const amrex::Real wind_speed_init,
    const amrex::Real yaw_init,
    const amrex::Real generator_efficiency)
{
    // Number of blades
    const auto& n_blades_wio = wio["assembly"]["number_of_blades"];
    if (n_blades_wio && (n_blades != n_blades_wio.as<int>())) {
        abort_run(
            "Number of blades in the YAML input to Kynema-FMB (" +
            std::to_string(n_blades_wio.as<int>()) + ") does not match " +
            "the expected number of blades from the Kynema-SGF input file "
            "(" +
            std::to_string(n_blades) + ").");
    }

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
            .SetSectionRefinement(section_refinement_blade);

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
                kynema_fmb::interfaces::components::ReferenceAxisOrientation::
                    Z);
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
                kynema_fmb::interfaces::components::ReferenceAxisOrientation::
                    Z);
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
        .SetSectionRefinement(section_refinement_tower)
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
            kynema_fmb::interfaces::components::ReferenceAxisOrientation::Z);
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
        const auto section = kynema_fmb::beams::GenerateHollowCircleSection(
            tower_diameter_grid[i], elastic_modulus, shear_modulus, density,
            tower_diameter_values[i], tower_wall_thickness_values[i],
            poisson_ratio);

        // Add section
        tower_builder.AddSection(
            section.position, section.M_star, section.C_star,
            kynema_fmb::interfaces::components::ReferenceAxisOrientation::Z);
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

std::array<int, 2> build_aero(
    kynema_fmb::interfaces::TurbineInterfaceBuilder& builder,
    const YAML::Node wio,
    const int n_blades,
    const bool do_tower_aero)
{
    //--------------------------------------------------------------------------
    // Build Aerodynamics
    //--------------------------------------------------------------------------

    auto airfoil_blade_map = std::vector(n_blades, 0UL);
    if (do_tower_aero) {
        airfoil_blade_map.emplace_back(1UL);
    }

    auto& aero_builder = builder.Aerodynamics()
                             .EnableAero()
                             .SetNumberOfAirfoils(do_tower_aero ? 2UL : 1UL)
                             .SetAirfoilToBladeMap(airfoil_blade_map);

    const auto& airfoil_io = wio["airfoils"];
    auto blade_aero_sections =
        std::vector<kynema_fmb::interfaces::components::AerodynamicSection>{};
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
        const auto aoa_copy = aoa;
        const auto cl = af["polars"][0]["re_sets"][0]["cl"]["values"]
                            .as<std::vector<amrex::Real>>();
        const auto cd = af["polars"][0]["re_sets"][0]["cd"]["values"]
                            .as<std::vector<amrex::Real>>();
        const auto cm = af["polars"][0]["re_sets"][0]["cm"]["values"]
                            .as<std::vector<amrex::Real>>();

        blade_aero_sections.emplace_back(
            id, static_cast<double>(s), static_cast<double>(chord),
            static_cast<double>(section_offset_x),
            static_cast<double>(section_offset_y),
            static_cast<double>(aerodynamic_center), static_cast<double>(twist),
            std::vector<double>(aoa_copy.begin(), aoa_copy.end()),
            std::vector<double>(cl.begin(), cl.end()),
            std::vector<double>(cd.begin(), cd.end()),
            std::vector<double>(cm.begin(), cm.end()));
        ++id;
    }

    aero_builder.SetAirfoilSections(0UL, blade_aero_sections);

    // Tower airfoil sections (simple circular cylinder)
    auto tower_aero_sections =
        std::vector<kynema_fmb::interfaces::components::AerodynamicSection>{};
    if (do_tower_aero) {
        const auto& tower_os = wio["components"]["tower"]["outer_shape"];
        const auto aoa = std::vector<amrex::Real>{-180.0_rt, 180.0_rt};
        const auto cl = std::vector<amrex::Real>{0.0_rt, 0.0_rt};
        const auto cm = std::vector<amrex::Real>{0.0_rt, 0.0_rt};
        const auto twist = 0.0_rt;
        const auto section_offset_x = 0.0_rt;
        const auto section_offset_y = 0.0_rt;
        const auto aerodynamic_center = 0.5_rt;
        auto cd_vec = tower_os["cd"]["values"].as<std::vector<amrex::Real>>();
        auto s_cd = tower_os["cd"]["grid"].as<std::vector<amrex::Real>>();
        if (tower_os.size() != tower_os["outer_diameter"]["values"]
                                   .as<std::vector<amrex::Real>>()
                                   .size()) {
            print_message(
                "WARNING: number of tower aero sections does not "
                "match number of cd tower data. Linearly "
                "interpolating cd to tower aero sections.\n");
        }
        id = 0UL;
        for (const auto& chord : tower_os["outer_diameter"]["values"]
                                     .as<std::vector<amrex::Real>>()) {
            const auto s = tower_os["outer_diameter"]["grid"]
                               .as<std::vector<amrex::Real>>()[id];
            const auto cd = interp_linear(
                s_cd.data(), cd_vec.data(), static_cast<int>(s_cd.size()), s);
            const auto cd_vec_wrapped = std::vector<amrex::Real>{cd, cd};
            tower_aero_sections.emplace_back(
                id, static_cast<double>(s), static_cast<double>(chord),
                static_cast<double>(section_offset_x),
                static_cast<double>(section_offset_y),
                static_cast<double>(aerodynamic_center),
                static_cast<double>(twist),
                std::vector<double>(aoa.begin(), aoa.end()),
                std::vector<double>(cl.begin(), cl.end()),
                std::vector<double>(
                    cd_vec_wrapped.begin(), cd_vec_wrapped.end()),
                std::vector<double>(cm.begin(), cm.end()));
            ++id;
        }
        aero_builder.SetAirfoilSections(1UL, tower_aero_sections);
    }

    return std::array<int, 2>{
        (int)blade_aero_sections.size(), (int)tower_aero_sections.size()};
}

void build_controller(
    kynema_fmb::interfaces::TurbineInterfaceBuilder& builder,
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
                                  .SetOutputFilePath(controller_output_file);
}
} // namespace

TurbineInterfacePtr create_turbine_interface(const TurbineSetup& setup)
{
    auto builder = kynema_fmb::interfaces::TurbineInterfaceBuilder{};

    builder.Solution()
        .EnableDynamicSolve()
        .SetTimeStep(setup.dt_ext)
        .SetDampingFactor(setup.damping_factor)
        .SetGravity({setup.gravity[0], setup.gravity[1], setup.gravity[2]})
        .SetMaximumNonlinearIterations(static_cast<size_t>(setup.nl_iter_max))
        .SetAbsoluteErrorTolerance(setup.abs_err_tol)
        .SetRelativeErrorTolerance(setup.rel_err_tol);

    const YAML::Node wio = YAML::LoadFile(setup.input_file);

    // Builds turbine, including blades, nacelle, and tower
    build_turbine(
        builder, wio, setup.num_blades, setup.num_blade_elem,
        setup.num_tower_elem, setup.section_refinement_blade,
        setup.section_refinement_tower, setup.rotational_speed,
        setup.generator_power, setup.wind_speed, setup.yaw,
        setup.generator_efficiency);

    auto n_aero_sections =
        build_aero(builder, wio, setup.num_blades, (setup.num_pts_tower != 0));

    if (n_aero_sections[0] != setup.num_pts_blade) {
        abort_run(
            "KynemaFMBIface: number of points per blade (from Kynema-SGF input "
            "file: " +
            std::to_string(setup.num_pts_blade) +
            ") does not match number of aerodynamic sections per blade (from "
            "Kynema-FMB input file: " +
            std::to_string(n_aero_sections[0]) + ").");
    }

    if (n_aero_sections[1] != setup.num_pts_tower) {
        abort_run(
            "KynemaFMBIface: number of tower points (from Kynema-SGF input "
            "file: " +
            std::to_string(setup.num_pts_tower) +
            ") does not match number of aerodynamic tower sections (from "
            "Kynema-FMB input file: " +
            std::to_string(n_aero_sections[1]) + ").");
    }

    if (setup.controller_input_file.size() > 0) {
        build_controller(
            builder, setup.controller_shared_lib_path,
            setup.controller_input_file, setup.label);
        // Should this enable controller checkpoint reading if it is a restart?
        // TODO? No Kynema-FMB examples to reference for this
    }

    // Create output
    builder.Outputs().SetOutputFilePath("kynema_fmb_" + setup.label);

    // Create checkpoint file
    if (!setup.checkpoint_file.empty()) {
        builder.Turbine().SetRestartFilePath(setup.checkpoint_file);
    }

    TurbineInterfacePtr ti(new TurbineInterface(
        builder.Solution().Input(), builder.Turbine().Input(),
        builder.Aerodynamics().Input(), builder.Controller().Input(),
        builder.Outputs().Config()));

    // Close file
    ti->CloseOutputFile();

    return ti;
}

void set_time_step_index(TurbineInterface& ti, const int time_index)
{
    ti.SetTimeStepIndex(time_index);
}

void get_hub_state(
    TurbineInterface& ti, float* hub_abs_pos, float* hub_rot_vel)
{
    const auto& hub_node = ti.Turbine().hub_node;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        hub_abs_pos[dir] = hub_node.position[dir];
        hub_rot_vel[dir] = hub_node.velocity[dir + AMREX_SPACEDIM];
    }
}

void write_checkpoint_file(TurbineInterface& ti, const std::string& fname)
{
    ti.WriteCheckpointFile(fname);
}

void calculate_motion(TurbineInterface& ti)
{
    ti.Aerodynamics().CalculateMotion(ti.GetHostState());
}

void calculate_aerodynamic_loads(TurbineInterface& ti, const amrex::Real density)
{
    ti.Aerodynamics().CalculateAerodynamicLoads(density);
}

void calculate_nodal_loads(TurbineInterface& ti)
{
    ti.Aerodynamics().CalculateNodalLoads();
}

bool step(TurbineInterface& ti) { return ti.Step(); }

void apply_controller(TurbineInterface& ti, const amrex::Real time)
{
    ti.ApplyController(time);
}

void write_output(TurbineInterface& ti)
{
    ti.OpenOutputFile();
    ti.WriteOutput();
    ti.CloseOutputFile();
}

} // namespace sgf_fmb
