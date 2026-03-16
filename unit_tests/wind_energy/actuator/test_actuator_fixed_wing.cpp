#include "aw_test_utils/MeshTest.H"
#include "test_act_utils.H"

#include "amr-wind/wind_energy/actuator/Actuator.H"
#include "amr-wind/wind_energy/actuator/ActuatorContainer.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "amr-wind/wind_energy/actuator/wing/fixed_wing_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {
namespace {
namespace act = amr_wind::actuator;

struct FixedWing : public act::FixedWing
{
    // cppcheck-suppress duplInheritedMember
    static std::string identifier() { return "TestFixedWing"; }
};

class ActFixedWingTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{32, 32, 32}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 16);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", 0.2_rt);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{-16.0_rt, -16.0_rt, -16.0_rt}};
            amrex::Vector<amrex::Real> probhi{{16.0_rt, 16.0_rt, 16.0_rt}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }

    void initialize_domain()
    {
        initialize_mesh();
        sim().repo().declare_field("actuator_src_term", 3, 0);
        auto& vel = sim().repo().declare_field("velocity", 3, 3);
        auto& density = sim().repo().declare_field("density", 1, 3);
        vel.setVal(10.0_rt, 0, 1, 3);
        density.setVal(1.0_rt);
        amr_wind::actuator::ActuatorContainer::ParticleType::NextID(1U);
    }

    void moving_wing_setup()
    {
        amrex::ParmParse pp_a("Actuator");
        pp_a.add("labels", (std::string) "F1");
        pp_a.add("type", (std::string) "TestFixedWingLine");
        amrex::ParmParse pp("Actuator.TestFixedWingLine");
        pp.add("num_points", 21);
        pp.add("epsilon", 1.0_rt);
        pp.add("pitch", 4.0_rt);
        pp.add("airfoil_table", m_afname);
        pp.add("airfoil_type", (std::string) "openfast");
        pp.add("motion_type", (std::string) "linear");
        pp.addarr(
            "velocity", amrex::Vector<amrex::Real>{1.0_rt, 0.5_rt, 0.7_rt});
        pp.addarr("span_locs", amrex::Vector<amrex::Real>{0.0_rt, 1.0_rt});
        pp.addarr("chord", amrex::Vector<amrex::Real>{2.0_rt, 2.0_rt});
        amrex::ParmParse pp_f("Actuator.F1");
        pp_f.addarr(
            "start", amrex::Vector<amrex::Real>{0.0_rt, -4.0_rt, 0.0_rt});
        pp_f.addarr("end", amrex::Vector<amrex::Real>{0.0_rt, 4.0_rt, 0.0_rt});
    }

    void pitching_wing_2d_setup()
    {
        amrex::ParmParse pp_a("Actuator");
        pp_a.add("labels", (std::string) "F1");
        pp_a.add("type", (std::string) "TestFixedWingLine");
        amrex::ParmParse pp("Actuator.TestFixedWingLine");
        pp.add("num_points", 2);
        pp.add("airfoil_table", m_afname);
        pp.add("airfoil_type", (std::string) "openfast");
        pp.addarr(
            "epsilon", amrex::Vector<amrex::Real>{1.0_rt, 2.0_rt, 1.0_rt});
        pp.addarr("span_locs", amrex::Vector<amrex::Real>{0.0_rt, 1.0_rt});
        pp.addarr("chord", amrex::Vector<amrex::Real>{2.0_rt, 2.0_rt});
        amrex::ParmParse pp_f("Actuator.F1");
        pp_f.add("pitch_timetable", m_ptname);
        pp_f.add("disable_spanwise_gaussian", true);
        pp_f.addarr(
            "start", amrex::Vector<amrex::Real>{0.0_rt, -4.0_rt, 0.0_rt});
        pp_f.addarr("end", amrex::Vector<amrex::Real>{0.0_rt, 4.0_rt, 0.0_rt});
    }

    const std::string m_afname = "airfoil.txt";
    const std::string m_ptname = "pitch.txt";
};

std::stringstream generate_openfast_airfoil()
{
    std::stringstream ss;
    ss << "!........................................ " << '\n';
    ss << "! Table of aerodynamics coefficients " << '\n';
    ss << "        6   NumAlf            ! Number of data lines in the "
          "following table "
       << '\n';
    ss << "!    Alpha      Cl      Cd        Cm  " << '\n';
    ss << "!    (deg)      (-)     (-)       (-) " << '\n';
    ss << "   -180.00    0.000   0.0407   0.0000 " << '\n';
    ss << "   -175.00    0.223   0.0507   0.0937 " << '\n';
    ss << "   -170.00    0.405   0.1055   0.1702 " << '\n';
    ss << "   -160.00    0.658   0.2982   0.2819 " << '\n';
    ss << "   -155.00    0.733   0.4121   0.3213 " << '\n';
    ss << "   -150.00    0.778   0.5308   0.3520 " << '\n';
    ss << "   -145.00    0.795   0.6503   0.3754 " << '\n';
    ss << "   -140.00    0.787   0.7672   0.3926 " << '\n';
    ss << "   -135.00    0.757   0.8785   0.4046 " << '\n';

    return ss;
}

void write_airfoil_file(const std::string& fname)
{
    std::ofstream os(fname);
    os << generate_openfast_airfoil().str();
}

void write_pitch_file(const std::string& fname)
{
    std::ofstream os(fname);
    // Write header
    os << "Time\tPitchAngle\n";
    // Write time table
    os << "0.0\t0.0\n";
    os << "0.1\t90.0\n";
}

} // namespace
} // namespace amr_wind_tests

namespace amr_wind::actuator {
namespace ops {

template <>
struct ReadInputsOp<::amr_wind_tests::FixedWing, ActSrcLine>
{
    void operator()(
        ::amr_wind_tests::FixedWing::DataType& data, const utils::ActParser& pp)
    {
        ReadInputsOp<::amr_wind::actuator::FixedWing, ActSrcLine> actual_op;
        EXPECT_NO_FATAL_FAILURE(actual_op(data, pp));

        const auto& meta = data.meta();
        // Do checks for each case
        if (meta.pitch_timetable_file.empty()) {
            // Linear motion case
            EXPECT_DOUBLE_EQ(meta.vel_tr.x(), 1.0_rt);
            EXPECT_DOUBLE_EQ(meta.vel_tr.y(), 0.5_rt);
            EXPECT_DOUBLE_EQ(meta.vel_tr.z(), 0.7_rt);
        } else {
            // Pitching wing 2D case
            EXPECT_DOUBLE_EQ(meta.pitch_table[0], 0.0_rt);
            EXPECT_DOUBLE_EQ(meta.pitch_table[1], 90.0_rt);
        }
    }
};

template <>
struct InitDataOp<::amr_wind_tests::FixedWing, ActSrcLine>
{
    void operator()(::amr_wind_tests::FixedWing::DataType& data)
    {
        InitDataOp<::amr_wind::actuator::FixedWing, ActSrcLine> actual_op;
        EXPECT_NO_FATAL_FAILURE(actual_op(data));
    }
};

template <>
struct ComputeForceOp<::amr_wind_tests::FixedWing, ActSrcLine>
{
    void operator()(::amr_wind_tests::FixedWing::DataType& data)
    {
        constexpr amrex::Real tol =
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e1_rt;
        const auto& meta = data.meta();
        const auto& grid = data.grid();
        ComputeForceOp<::amr_wind::actuator::FixedWing, ActSrcLine> actual_op;
        EXPECT_NO_FATAL_FAILURE(actual_op(data));
        const auto time = data.sim().time().new_time();
        // Do checks for each case
        if (meta.pitch_timetable_file.empty()) {
            // Check position (using mid point)
            if (!(time > 0.0_rt)) {
                EXPECT_NEAR(grid.pos[10].x(), 1.0_rt * 0.0_rt, tol);
                EXPECT_NEAR(grid.pos[10].y(), 0.5_rt * 0.0_rt, tol);
                EXPECT_NEAR(grid.pos[10].z(), 0.7_rt * 0.0_rt, tol);
            } else {
                EXPECT_NEAR(grid.pos[10].x(), 1.0_rt * 0.2_rt, tol);
                EXPECT_NEAR(grid.pos[10].y(), 0.5_rt * 0.2_rt, tol);
                EXPECT_NEAR(grid.pos[10].z(), 0.7_rt * 0.2_rt, tol);
            }
        } else {
            // Check pitch
            if (!(time > 0.0_rt)) {
                EXPECT_NEAR(meta.pitch, 0.0_rt, tol);
            } else {
                EXPECT_NEAR(meta.pitch, 90.0_rt, tol);
            }
            // const auto dummy = grid.orientation;
            // Check rotation of epsilon
            amrex::Real angle;
            if (!(time > 0.0_rt)) {
                angle = 0.0_rt;
            } else {
                angle = 90.0_rt;
            }
            const auto ref_mat = vs::quaternion(meta.end - meta.start, angle);
            for (int n = 0; n < 9; ++n) {
                EXPECT_NEAR(grid.orientation[0][n], ref_mat[n], tol);
                EXPECT_NEAR(grid.orientation[1][n], ref_mat[n], tol);
            }
            // Check 2D nature of force
            EXPECT_NEAR(grid.dcoord_flags[1], 0.0_rt, tol);
        }
    }
};

template <>
struct ProcessOutputsOp<::amr_wind_tests::FixedWing, ActSrcLine>
{
    ProcessOutputsOp(::amr_wind_tests::FixedWing::DataType& /**/) {}
    void operator()(::amr_wind_tests::FixedWing::DataType& /**/) {}
    void read_io_options(const utils::ActParser& /**/) {}
    void prepare_outputs(const std::string& /**/) {}
    void write_outputs() {}
};

} // namespace ops
template class ::amr_wind::actuator::
    ActModel<::amr_wind_tests::FixedWing, ActSrcLine>;
} // namespace amr_wind::actuator

namespace amr_wind_tests {

class ActPhysicsTest : public ::amr_wind::actuator::Actuator
{
public:
    explicit ActPhysicsTest(::amr_wind::CFDSim& sim)
        : ::amr_wind::actuator::Actuator(sim)
    {}

protected:
    void prepare_outputs() final {}
};

TEST_F(ActFixedWingTest, linear_motion) // pitch_table, 2D
{
    write_airfoil_file(m_afname);
    initialize_domain();
    moving_wing_setup();
    ActPhysicsTest act(sim());
    act.pre_init_actions();
    act.post_init_actions();
    act.pre_advance_work();
    sim().time().new_timestep();
    act.pre_advance_work();
    // Delete airfoil
    const char* fname = m_afname.c_str();
    {
        std::ifstream f(fname);
        if (f.good()) {
            remove(fname);
        }
        // Check that file is removed
        std::ifstream ff(fname);
        EXPECT_FALSE(ff.good());
    }
}

TEST_F(ActFixedWingTest, pitch_table_2D)
{
    write_airfoil_file(m_afname);
    write_pitch_file(m_ptname);
    initialize_domain();
    pitching_wing_2d_setup();
    ActPhysicsTest act(sim());
    act.pre_init_actions();
    act.post_init_actions();
    act.pre_advance_work();
    sim().time().new_timestep();
    act.pre_advance_work();
    // Delete airfoil file
    {
        const char* fname = m_afname.c_str();
        std::ifstream f(fname);
        if (f.good()) {
            remove(fname);
        }
        // Check that file is removed
        std::ifstream ff(fname);
        EXPECT_FALSE(ff.good());
    }
    // Delete pitch file
    {
        const char* fname = m_ptname.c_str();
        std::ifstream f(fname);
        if (f.good()) {
            remove(fname);
        }
        // Check that file is removed
        std::ifstream ff(fname);
        EXPECT_FALSE(ff.good());
    }
}

} // namespace amr_wind_tests
