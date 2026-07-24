#include "src/wind_energy/actuator/drone/Drone.H"
#include "src/wind_energy/actuator/drone/drone_ops.H"
#include "src/wind_energy/actuator/Actuator.H"
#include "src/wind_energy/actuator/ActuatorContainer.H"
#include "src/wind_energy/actuator/ActuatorModel.H"
#include "src/utilities/constants.H"
#include "src/utilities/ncutils/nc_interface.H"
#include "ks_test_utils/MeshTest.H"

#include "gtest/gtest.h"

#include <cmath>
#include <filesystem>
#include <fstream>

#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {
namespace {

using kynema_sgf::actuator::drone::rotor_body_offsets;
using kynema_sgf::actuator::drone::uniform_arm_angles;
constexpr amrex::Real test_tol = kynema_sgf::constants::TIGHT_TOL;

class DroneActuatorTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();
        amrex::ParmParse pp_amr("amr");
        pp_amr.add("max_level", 0);
        pp_amr.add("max_grid_size", 16);
        pp_amr.addarr("n_cell", amrex::Vector<int>{16, 16, 16});

        amrex::ParmParse pp_time("time");
        pp_time.add("fixed_dt", 1.0e-4_rt);

        amrex::ParmParse pp_geom("geometry");
        pp_geom.addarr(
            "prob_lo", amrex::Vector<amrex::Real>{-0.2_rt, -0.2_rt, -0.2_rt});
        pp_geom.addarr(
            "prob_hi", amrex::Vector<amrex::Real>{0.2_rt, 0.2_rt, 0.2_rt});
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

    void populate_inputs()
    {
        amrex::ParmParse pp_a("Actuator");
        pp_a.add("labels", std::string("D1"));

        amrex::ParmParse pp_d("Actuator.Drone");
        pp_d.add("num_rotors", 4);
        pp_d.add("arm_length", 0.075_rt);
        pp_d.add("arm_phase_degrees", 45.0_rt);
        pp_d.add("rotor_diameter", 0.05_rt);
        pp_d.add("num_blades", 2);
        pp_d.add("root_radius_fraction", 0.18_rt);
        pp_d.add("epsilon_chord", 0.5_rt);
        pp_d.add("airfoil_table", m_airfoil_file);
        pp_d.add("airfoil_type", std::string("openfast"));
        pp_d.addarr("span_locs", amrex::Vector<amrex::Real>{0.0_rt, 1.0_rt});
        pp_d.addarr("chord", amrex::Vector<amrex::Real>{0.01_rt, 0.006_rt});
        pp_d.addarr("twist", amrex::Vector<amrex::Real>{12.0_rt, 4.0_rt});

        amrex::ParmParse pp_i("Actuator.D1");
        pp_i.add("type", std::string("Drone"));
        pp_i.addarr(
            "arm_length",
            amrex::Vector<amrex::Real>{0.075_rt, 0.075_rt, 0.075_rt, 0.075_rt});
        pp_i.addarr(
            "arm_angles_degrees",
            amrex::Vector<amrex::Real>{0.0_rt, 90.0_rt, 180.0_rt, 270.0_rt});
        pp_i.addarr(
            "center", amrex::Vector<amrex::Real>{0.0_rt, 0.0_rt, 0.0_rt});
        pp_i.addarr(
            "translation_velocity",
            amrex::Vector<amrex::Real>{0.1_rt, -0.2_rt, 0.3_rt});
        pp_i.addarr(
            "rotor_omegas", amrex::Vector<amrex::Real>{
                                2500.0_rt, -2500.0_rt, 2500.0_rt, -2500.0_rt});
        pp_i.addarr(
            "mirror_blades",
            amrex::Vector<std::string>{"false", "true", "false", "true"});

        amrex::ParmParse pp_s("Actuator.ActuatorSector");
        pp_s.add("rotor_diameter", 0.05_rt);
        pp_s.add("num_blades", 2);
        pp_s.add("root_radius_fraction", 0.18_rt);
        pp_s.add("epsilon_chord", 0.5_rt);
        pp_s.add("airfoil_table", m_airfoil_file);
        pp_s.add("airfoil_type", std::string("openfast"));
        pp_s.addarr("span_locs", amrex::Vector<amrex::Real>{0.0_rt, 1.0_rt});
        pp_s.addarr("chord", amrex::Vector<amrex::Real>{0.01_rt, 0.006_rt});
        pp_s.addarr("twist", amrex::Vector<amrex::Real>{12.0_rt, 4.0_rt});

        const amrex::Real offset = 0.075_rt / std::sqrt(2.0_rt);
        const amrex::Vector<amrex::Vector<amrex::Real>> centers{
            {offset, offset, 0.0_rt},
            {-offset, offset, 0.0_rt},
            {-offset, -offset, 0.0_rt},
            {offset, -offset, 0.0_rt}};
        for (int i = 0; i < 4; ++i) {
            amrex::ParmParse pp_r("Actuator.R" + std::to_string(i + 1));
            pp_r.add("omega", (i % 2 == 0) ? 2500.0_rt : -2500.0_rt);
            if (i % 2 != 0) {
                pp_r.addarr(
                    "twist", amrex::Vector<amrex::Real>{-12.0_rt, -4.0_rt});
            }
            pp_r.addarr("center", centers[i]);
            pp_r.addarr(
                "translation_velocity",
                amrex::Vector<amrex::Real>{0.1_rt, -0.2_rt, 0.3_rt});
            pp_r.addarr(
                "rotor_normal",
                amrex::Vector<amrex::Real>{0.0_rt, 0.0_rt, 1.0_rt});
        }
    }

    void write_airfoil() const
    {
        std::ofstream os(m_airfoil_file);
        os << "! test polar\n5 NumAlf\n! Alpha Cl Cd Cm\n! deg - - -\n"
              "-180 0 0.04 0\n-10 -0.5 0.02 0\n0 0 0.01 0\n"
              "10 0.5 0.02 0\n180 0 0.04 0\n";
    }

    const std::string m_airfoil_file{"drone_airfoil.txt"};
};

class DronePhysicsTest : public kynema_sgf::actuator::Actuator
{
public:
    explicit DronePhysicsTest(kynema_sgf::CFDSim& sim) : Actuator(sim) {}

protected:
    void prepare_outputs() override {}
};

TEST(DroneGeometry, plus_layout)
{
    const auto angles = uniform_arm_angles(4, 0.0_rt);
    const auto offsets =
        rotor_body_offsets({2.0_rt, 2.0_rt, 2.0_rt, 2.0_rt}, angles);

    ASSERT_EQ(offsets.size(), 4U);
    EXPECT_NEAR(offsets[0].x(), 2.0_rt, test_tol);
    EXPECT_NEAR(offsets[0].y(), 0.0_rt, test_tol);
    EXPECT_NEAR(offsets[1].x(), 0.0_rt, test_tol);
    EXPECT_NEAR(offsets[1].y(), 2.0_rt, test_tol);
    EXPECT_NEAR(offsets[2].x(), -2.0_rt, test_tol);
    EXPECT_NEAR(offsets[3].y(), -2.0_rt, test_tol);
}

TEST(DroneGeometry, unequal_irregular_arms)
{
    const kynema_sgf::actuator::RealList lengths{1.0_rt, 2.0_rt, 3.0_rt};
    const kynema_sgf::actuator::RealList angles{0.0_rt, 90.0_rt, 225.0_rt};
    const auto offsets = rotor_body_offsets(lengths, angles);

    EXPECT_NEAR(offsets[0].x(), 1.0_rt, test_tol);
    EXPECT_NEAR(offsets[1].y(), 2.0_rt, test_tol);
    EXPECT_NEAR(offsets[2].x(), -3.0_rt / std::sqrt(2.0_rt), test_tol);
    EXPECT_NEAR(offsets[2].y(), -3.0_rt / std::sqrt(2.0_rt), test_tol);
}

TEST(DroneGeometry, body_orientation_defaults_to_identity)
{
    const auto rotation = kynema_sgf::actuator::drone::body_rotation(
        kynema_sgf::vs::Vector::zero());
    const auto value =
        rotation & kynema_sgf::vs::Vector{1.0_rt, 2.0_rt, 3.0_rt};
    EXPECT_NEAR(value.x(), 1.0_rt, test_tol);
    EXPECT_NEAR(value.y(), 2.0_rt, test_tol);
    EXPECT_NEAR(value.z(), 3.0_rt, test_tol);
}

TEST_F(DroneActuatorTest, composite_lifecycle)
{
    write_airfoil();
    initialize_domain();
    populate_inputs();

    DronePhysicsTest actuator(sim());
    actuator.pre_init_actions();
    auto* drone = dynamic_cast<kynema_sgf::actuator::ActModel<
        kynema_sgf::actuator::Drone, kynema_sgf::actuator::ActSrcDrone>*>(
        &actuator.get_act(0));
    ASSERT_NE(drone, nullptr);
    ASSERT_EQ(drone->meta().rotors.size(), 4U);
    EXPECT_NEAR(
        drone->meta().rotors[0]->data.meta().center.x(),
        0.075_rt / std::sqrt(2.0_rt), test_tol);
    EXPECT_NEAR(
        drone->meta().rotors[0]->data.meta().center.y(),
        0.075_rt / std::sqrt(2.0_rt), test_tol);

    actuator.post_init_actions();
    actuator.pre_advance_work();

    remove(m_airfoil_file.c_str());
}

TEST_F(DroneActuatorTest, matches_equivalent_standalone_sectors)
{
    namespace act = kynema_sgf::actuator;
    constexpr amrex::Real tol = kynema_sgf::constants::TIGHT_TOL;

    write_airfoil();
    initialize_domain();
    populate_inputs();

    act::ActModel<act::Drone, act::ActSrcDrone> drone(sim(), "D1", 0);
    drone.read_inputs(act::utils::ActParser("Actuator.Drone", "Actuator.D1"));
    drone.init_actuator_source();
    ASSERT_EQ(drone.meta().rotors.size(), 4U);

    for (int i = 0; i < 4; ++i) {
        const std::string label = "R" + std::to_string(i + 1);
        act::ActModel<act::ActuatorSector, act::ActSrcSector> sector(
            sim(), label, i);
        sector.read_inputs(
            act::utils::ActParser(
                "Actuator.ActuatorSector", "Actuator." + label));
        sector.init_actuator_source();

        const auto& drone_data = drone.meta().rotors[i]->data;
        const auto& drone_meta = drone_data.meta();
        const auto& sector_meta = sector.meta();
        const auto& drone_grid = drone_data.grid();
        const auto& sector_grid = sector.grid();

        EXPECT_EQ(drone_meta.num_blades, sector_meta.num_blades);
        EXPECT_NEAR(drone_meta.rotor_diameter, sector_meta.rotor_diameter, tol);
        EXPECT_NEAR(drone_meta.rotor_radius, sector_meta.rotor_radius, tol);
        EXPECT_NEAR(drone_meta.root_radius, sector_meta.root_radius, tol);
        EXPECT_NEAR(drone_meta.omega, sector_meta.omega, tol);

        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            EXPECT_NEAR(drone_meta.center[n], sector_meta.center[n], tol);
            EXPECT_NEAR(
                drone_meta.translation_velocity[n],
                sector_meta.translation_velocity[n], tol);
            EXPECT_NEAR(
                drone_meta.rotor_normal[n], sector_meta.rotor_normal[n], tol);
        }

        ASSERT_EQ(drone_meta.radius.size(), sector_meta.radius.size());
        EXPECT_EQ(drone_meta.dr.size(), sector_meta.dr.size());
        EXPECT_EQ(drone_meta.chord.size(), sector_meta.chord.size());
        EXPECT_EQ(drone_meta.twist.size(), sector_meta.twist.size());
        EXPECT_EQ(
            drone_meta.epsilon_profile.size(),
            sector_meta.epsilon_profile.size());
        EXPECT_EQ(drone_meta.theta_counts, sector_meta.theta_counts);
        for (int j = 0; j < static_cast<int>(drone_meta.radius.size()); ++j) {
            EXPECT_NEAR(drone_meta.radius[j], sector_meta.radius[j], tol);
            EXPECT_NEAR(drone_meta.dr[j], sector_meta.dr[j], tol);
            EXPECT_NEAR(drone_meta.chord[j], sector_meta.chord[j], tol);
            EXPECT_NEAR(drone_meta.twist[j], sector_meta.twist[j], tol);
            EXPECT_NEAR(
                drone_meta.epsilon_profile[j], sector_meta.epsilon_profile[j],
                tol);
        }

        EXPECT_EQ(drone_grid.pos.size(), sector_grid.pos.size());
        EXPECT_EQ(drone_grid.force.size(), sector_grid.force.size());
        EXPECT_EQ(drone_grid.epsilon.size(), sector_grid.epsilon.size());
        ASSERT_EQ(drone_grid.vel_pos.size(), sector_grid.vel_pos.size());
        EXPECT_EQ(drone_grid.vel.size(), sector_grid.vel.size());
        EXPECT_EQ(drone_grid.density.size(), sector_grid.density.size());
        for (int j = 0; j < static_cast<int>(drone_grid.vel_pos.size()); ++j) {
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                EXPECT_NEAR(
                    drone_grid.vel_pos[j][n], sector_grid.vel_pos[j][n], tol);
            }
        }
    }

    remove(m_airfoil_file.c_str());
}

TEST_F(
    DroneActuatorTest, counter_rotating_rotors_have_same_axial_force_direction)
{
    namespace act = kynema_sgf::actuator;

    write_airfoil();
    initialize_domain();
    populate_inputs();

    act::ActModel<act::Drone, act::ActSrcDrone> drone(sim(), "D1", 0);
    drone.read_inputs(act::utils::ActParser("Actuator.Drone", "Actuator.D1"));
    drone.init_actuator_source();
    ASSERT_EQ(drone.meta().rotors.size(), 4U);

    amrex::Real reference_force_normal = 0.0_rt;
    for (int i = 0; i < 4; ++i) {
        auto& rotor_data = drone.meta().rotors[i]->data;
        act::ops::ComputeForceOp<act::ActuatorSector, act::ActSrcSector>()(
            rotor_data);
        const auto& rotor_meta = rotor_data.meta();
        const amrex::Real force_normal =
            rotor_meta.integrated_force & rotor_meta.rotor_normal;
        ASSERT_NE(force_normal, 0.0_rt);
        if (i == 0) {
            reference_force_normal = force_normal;
        } else {
            EXPECT_GT(force_normal * reference_force_normal, 0.0_rt);
        }
    }

    remove(m_airfoil_file.c_str());
}

#ifdef KYNEMA_SGF_USE_NETCDF
TEST_F(DroneActuatorTest, writes_rotors_as_nested_netcdf_groups)
{
    namespace act = kynema_sgf::actuator;

    write_airfoil();
    initialize_domain();
    populate_inputs();

    act::ActModel<act::Drone, act::ActSrcDrone> drone(sim(), "D1", 0);
    drone.read_inputs(act::utils::ActParser("Actuator.Drone", "Actuator.D1"));
    drone.init_actuator_source();
    drone.prepare_outputs(".");
    drone.write_outputs();

    auto ncf = ncutils::NCFile::open("D1.nc");
    EXPECT_EQ(ncf.dim("num_time_steps").len(), 1U);
    ASSERT_TRUE(ncf.has_group("D1"));
    auto drone_group = ncf.group("D1");
    EXPECT_TRUE(drone_group.has_var("force"));
    EXPECT_TRUE(drone_group.has_var("moment"));
    for (int i = 0; i < 4; ++i) {
        const std::string rotor_name = "R" + std::to_string(i + 1);
        ASSERT_TRUE(drone_group.has_group(rotor_name));
        auto rotor_group = drone_group.group(rotor_name);
        EXPECT_TRUE(rotor_group.has_var("thrust"));
        EXPECT_TRUE(rotor_group.has_var("torque"));
        EXPECT_TRUE(rotor_group.has_var("blade_force"));
        EXPECT_TRUE(rotor_group.has_var("aoa"));
        EXPECT_EQ(rotor_group.var("time").shape().front(), 1U);
    }
    ncf.close();

    EXPECT_FALSE(std::filesystem::exists("D1.R1.nc"));
    remove("D1.nc");
    remove(m_airfoil_file.c_str());
}
#endif

} // namespace
} // namespace kynema_sgf_tests
