#include <limits>
#include <string>
#include "ks_test_utils/AmrexTest.H"
#include "src/actuator/FLLC.H"
#include "src/actuator/FLLCOp.H"
#include "src/utilities/linear_interpolation.H"
#include "AMReX_Exception.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::actuator {

TEST(TestFLLCData, data_initializes_with_cviews)
{
    const int num_points = 10;

    VecList dummy_vec(num_points);
    RealList dummy_real(num_points);
    TensorList dummy_tensor(num_points);

    FLLCData data;

    ComponentView view;
    view.pos = ::kynema_sgf::utils::slice(dummy_vec, 0, num_points);
    view.force = ::kynema_sgf::utils::slice(dummy_vec, 0, num_points);
    view.epsilon = ::kynema_sgf::utils::slice(dummy_vec, 0, num_points);
    view.orientation = ::kynema_sgf::utils::slice(dummy_tensor, 0, num_points);
    view.chord = ::kynema_sgf::utils::slice(dummy_real, 0, num_points);
    view.vel_rel = ::kynema_sgf::utils::slice(dummy_vec, 0, num_points);

    data.nonuniform = false;

    for (int ip = 0; ip < num_points; ++ip) {
        view.pos[ip] = {0, static_cast<amrex::Real>(ip), 0};
        view.chord[ip] = 1.;
    }

    fllc_init(data, view, 1.0_rt);

    ASSERT_EQ(num_points, data.les_velocity.size());
    ASSERT_EQ(num_points, data.optimal_velocity.size());
    ASSERT_EQ(num_points, data.lift.size());
    ASSERT_EQ(num_points, data.grad_lift.size());
}

TEST(TestFLLCData, data_initializes_with_cviews_nonuniform)
{
    const int num_points = 5;

    VecList dummy_vec(num_points);
    RealList dummy_real(num_points);
    TensorList dummy_tensor(num_points);

    FLLCData data;

    ComponentView view;
    view.pos = ::kynema_sgf::utils::slice(dummy_vec, 0, num_points);
    view.force = ::kynema_sgf::utils::slice(dummy_vec, 0, num_points);
    view.epsilon = ::kynema_sgf::utils::slice(dummy_vec, 0, num_points);
    view.orientation = ::kynema_sgf::utils::slice(dummy_tensor, 0, num_points);
    view.chord = ::kynema_sgf::utils::slice(dummy_real, 0, num_points);
    view.vel_rel = ::kynema_sgf::utils::slice(dummy_vec, 0, num_points);

    data.nonuniform = true;
    data.eps_dr = 1.1_rt;

    for (int ip = 0; ip < num_points; ++ip) {
        view.pos[ip] = {0, static_cast<amrex::Real>(ip), 0};
        view.chord[ip] = 1.;
    }

    fllc_init(data, view, 1.0_rt);

    const int npts_r = static_cast<int>(data.nonuniform_r.size());
    const int npts_dr = static_cast<int>(data.nonuniform_dr.size());

    ASSERT_EQ(npts_r, npts_dr);
}

namespace {

/** Straight blade of a few points for the tests of the correction
 *
 *  Force points one unit apart along y, starting at start, with a unit chord
 *  and a relative velocity along x. The velocity points are given by the
 *  test: they need not be the force points.
 */
struct FLLCTestBlade
{
    FLLCTestBlade(const int npts, const int npts_vel, const vs::Vector& start)
        : pos(npts)
        , force(npts, vs::Vector::zero())
        , epsilon(npts)
        , vel_rel(npts, vs::Vector{1.0_rt, 0.0_rt, 0.0_rt})
        , orientation(npts)
        , chord(npts, 1.0_rt)
        , vel_pos(npts_vel)
        , vel(npts_vel, vs::Vector::zero())
    {
        for (int ip = 0; ip < npts; ++ip) {
            const vs::Vector along{
                0.0_rt, static_cast<amrex::Real>(ip), 0.0_rt};
            pos[ip] = start + along;
        }
        view.pos = ::kynema_sgf::utils::slice(pos, 0, npts);
        view.force = ::kynema_sgf::utils::slice(force, 0, npts);
        view.epsilon = ::kynema_sgf::utils::slice(epsilon, 0, npts);
        view.orientation = ::kynema_sgf::utils::slice(orientation, 0, npts);
        view.chord = ::kynema_sgf::utils::slice(chord, 0, npts);
        view.vel_rel = ::kynema_sgf::utils::slice(vel_rel, 0, npts);
        view.vel_pos = ::kynema_sgf::utils::slice(vel_pos, 0, npts_vel);
        view.vel = ::kynema_sgf::utils::slice(vel, 0, npts_vel);
    }

    //! Put velocity point ip at the given distances along the blade
    void set_vel_span(const RealList& span, const vs::Vector& shift)
    {
        for (int ip = 0; ip < static_cast<int>(vel_pos.size()); ++ip) {
            const vs::Vector along{0.0_rt, span[ip], 0.0_rt};
            vel_pos[ip] = pos[0] + shift + along;
        }
    }

    VecList pos;
    VecList force;
    VecList epsilon;
    VecList vel_rel;
    TensorList orientation;
    RealList chord;
    VecList vel_pos;
    VecList vel;
    ComponentView view;
};

} // namespace

TEST(TestFLLCOp, kernel_keeps_precision_at_small_distances)
{
    const int num_points = 5;
    const RealList span{0.0_rt, 1.0_rt, 2.0_rt, 3.0_rt, 4.0_rt};

    // The same blade twice. Epsilon over the spacing of the non-uniform
    // points is 1, so these points are the force points.
    FLLCTestBlade blade_ref(num_points, num_points, vs::Vector::zero());
    FLLCTestBlade blade_shift(num_points, num_points, vs::Vector::zero());
    blade_ref.set_vel_span(span, vs::Vector::zero());
    blade_shift.set_vel_span(span, vs::Vector::zero());
    FLLCData data_ref;
    FLLCData data_shift;
    for (auto* data : {&data_ref, &data_shift}) {
        data->epsilon = 2.0_rt;
        data->relaxation_factor = 1.0_rt;
        data->nonuniform = true;
        data->eps_dr = 1.0_rt;
    }
    for (int ip = 0; ip < num_points; ++ip) {
        blade_ref.force[ip] = vs::Vector{0.0_rt, 0.0_rt, 1.0_rt};
        blade_shift.force[ip] = vs::Vector{0.0_rt, 0.0_rt, 1.0_rt};
    }
    fllc_init(data_ref, blade_ref.view, 1.0_rt);
    fllc_init(data_shift, blade_shift.view, 1.0_rt);
    ASSERT_EQ(static_cast<int>(data_ref.nonuniform_r.size()), num_points);

    // Move one non-uniform point off its force point by much less than
    // epsilon, as the round-off of the point builder does for a ratio of 3
    data_shift.nonuniform_r[2] += 1.0e-9_rt;

    FLLCOp()(blade_ref.view, data_ref);
    FLLCOp()(blade_shift.view, data_shift);

    // The correction cannot depend on a shift of 1e-9: the kernel is
    // continuous at a distance of 0
    for (int ip = 0; ip < num_points; ++ip) {
        const auto du_ref = data_ref.correction_velocity[ip].z();
        const auto du_shift = data_shift.correction_velocity[ip].z();
        EXPECT_GT(std::abs(du_ref), 1.0e-3_rt);
        EXPECT_NEAR(du_shift, du_ref, 1.0e-5_rt * std::abs(du_ref));
    }
}

TEST(TestFLLCOp, no_correction_when_force_is_spread_with_optimal_epsilon)
{
    const int num_points = 5;
    const RealList span{0.0_rt, 1.0_rt, 2.0_rt, 3.0_rt, 4.0_rt};

    // The optimal epsilon (chord times epsilon_chord = 1) is larger than
    // epsilon (0.5): the actuator spreads the force with the optimal one,
    // and there is nothing to correct
    for (const auto type :
         {FLLCType::VariableChord, FLLCType::ConstantChord}) {
        FLLCTestBlade blade(num_points, num_points, vs::Vector::zero());
        blade.set_vel_span(span, vs::Vector::zero());
        FLLCData data;
        data.correction_type = type;
        data.epsilon = 0.5_rt;
        data.relaxation_factor = 1.0_rt;
        data.nonuniform = true;
        data.eps_dr = 1.0_rt;
        for (int ip = 0; ip < num_points; ++ip) {
            blade.force[ip] = vs::Vector{0.0_rt, 0.0_rt, 1.0_rt};
            blade.vel[ip] = vs::Vector{8.0_rt, 0.0_rt, 0.0_rt};
        }
        fllc_init(data, blade.view, 1.0_rt);

        FLLCOp()(blade.view, data);

        for (int ip = 0; ip < num_points; ++ip) {
            EXPECT_EQ(vs::mag(data.correction_velocity[ip]), 0.0_rt);
            EXPECT_EQ(blade.vel[ip].x(), 8.0_rt);
            EXPECT_EQ(blade.vel[ip].z(), 0.0_rt);
        }
    }
}

namespace {

//! Force and settings shared by the tests of the mapping between the points
void init_mapping_test(FLLCTestBlade& blade, FLLCData& data)
{
    const RealList force_z{1.0_rt, 2.0_rt, 4.0_rt, 3.0_rt, 1.0_rt};
    for (int ip = 0; ip < static_cast<int>(blade.force.size()); ++ip) {
        blade.force[ip] = vs::Vector{0.0_rt, 0.0_rt, force_z[ip]};
    }
    data.epsilon = 2.0_rt;
    data.relaxation_factor = 1.0_rt;
    data.nonuniform = true;
    data.eps_dr = 1.0_rt;
    fllc_init(data, blade.view, 1.0_rt);
}

//! Correction at a distance along the blade, from the force points
amrex::Real
correction_at(const FLLCData& data, const int npts, const amrex::Real span)
{
    RealList force_span(npts);
    RealList du_z(npts);
    for (int ip = 0; ip < npts; ++ip) {
        force_span[ip] = static_cast<amrex::Real>(ip);
        du_z[ip] = data.correction_velocity[ip].z();
    }
    return ::kynema_sgf::interp::linear(force_span, du_z, span);
}

} // namespace

TEST(TestFLLCOp, same_points_at_different_instants)
{
    const int num_points = 5;
    const RealList span{0.0_rt, 1.0_rt, 2.0_rt, 3.0_rt, 4.0_rt};

    // The velocity points are the force points at another instant of a
    // moving blade: the same distances along the blade, another place
    FLLCTestBlade blade(num_points, num_points, vs::Vector::zero());
    blade.set_vel_span(span, vs::Vector{0.3_rt, 0.0_rt, 0.2_rt});
    FLLCData data;
    init_mapping_test(blade, data);

    FLLCOp()(blade.view, data);

    // Velocity point k takes the correction of force point k
    for (int ip = 0; ip < num_points; ++ip) {
        const auto du = data.correction_velocity[ip].z();
        EXPECT_GT(std::abs(du), 1.0e-3_rt);
        EXPECT_EQ(blade.vel[ip].z(), du);
        EXPECT_EQ(blade.vel[ip].x(), 0.0_rt);
    }
}

TEST(TestFLLCOp, same_number_of_points_at_other_positions)
{
    const int num_points = 5;
    const RealList span_vel{0.0_rt, 0.6_rt, 1.7_rt, 3.2_rt, 4.0_rt};
    const amrex::Real tol =
        100.0_rt * std::numeric_limits<amrex::Real>::epsilon();

    // As many velocity points as force points, but not at the force points
    FLLCTestBlade blade(num_points, num_points, vs::Vector::zero());
    blade.set_vel_span(span_vel, vs::Vector::zero());
    FLLCData data;
    init_mapping_test(blade, data);

    FLLCOp()(blade.view, data);

    // Every velocity point takes the correction at its own distance along
    // the blade
    for (int ip = 0; ip < num_points; ++ip) {
        const auto gold = correction_at(data, num_points, span_vel[ip]);
        EXPECT_NEAR(blade.vel[ip].z(), gold, tol) << "velocity point " << ip;
    }
}

TEST(TestFLLCOp, no_force_keeps_velocity_at_mismatched_points)
{
    const int num_points = 5;
    const int num_vel_points = 3;
    const RealList span_vel{0.0_rt, 1.5_rt, 4.0_rt};
    const RealList vel_x{1.0_rt, 5.0_rt, 2.0_rt};

    // A blade away from the origin with fewer velocity points than force
    // points, and no force
    FLLCTestBlade blade(
        num_points, num_vel_points, vs::Vector{10.0_rt, 5.0_rt, 0.0_rt});
    blade.set_vel_span(span_vel, vs::Vector::zero());
    for (int ip = 0; ip < num_vel_points; ++ip) {
        blade.vel[ip] = vs::Vector{vel_x[ip], 0.0_rt, 0.0_rt};
    }
    FLLCData data;
    data.epsilon = 2.0_rt;
    data.relaxation_factor = 1.0_rt;
    fllc_init(data, blade.view, 1.0_rt);

    FLLCOp()(blade.view, data);

    // No force, no correction: the sampled velocity is not to be touched
    for (int ip = 0; ip < num_vel_points; ++ip) {
        EXPECT_EQ(blade.vel[ip].x(), vel_x[ip]) << "velocity point " << ip;
        EXPECT_EQ(blade.vel[ip].y(), 0.0_rt);
        EXPECT_EQ(blade.vel[ip].z(), 0.0_rt);
    }
}

TEST(TestFLLCOp, correction_follows_span_at_mismatched_points)
{
    const int num_points = 5;
    const int num_vel_points = 3;
    const RealList span_vel{0.0_rt, 1.5_rt, 4.0_rt};
    const RealList vel_x{1.0_rt, 5.0_rt, 2.0_rt};
    const amrex::Real tol =
        100.0_rt * std::numeric_limits<amrex::Real>::epsilon();

    // A blade away from the origin with fewer velocity points than force
    // points
    FLLCTestBlade blade(
        num_points, num_vel_points, vs::Vector{10.0_rt, 5.0_rt, 0.0_rt});
    blade.set_vel_span(span_vel, vs::Vector::zero());
    for (int ip = 0; ip < num_vel_points; ++ip) {
        blade.vel[ip] = vs::Vector{vel_x[ip], 0.0_rt, 0.0_rt};
    }
    FLLCData data;
    init_mapping_test(blade, data);

    FLLCOp()(blade.view, data);

    // The correction is normal to the relative velocity (along x): the
    // sampled x velocity stays, and z takes the correction at the distance
    // of the velocity point along the blade
    for (int ip = 0; ip < num_vel_points; ++ip) {
        const auto gold = correction_at(data, num_points, span_vel[ip]);
        EXPECT_EQ(blade.vel[ip].x(), vel_x[ip]) << "velocity point " << ip;
        EXPECT_NEAR(blade.vel[ip].z(), gold, tol) << "velocity point " << ip;
    }
}

namespace {

//! Inputs of an actuator with the correction on, for the tests of fllc_parse
class FLLCParseTest : public ::kynema_sgf_tests::AmrexTest
{
protected:
    void SetUp() override
    {
        ::kynema_sgf_tests::AmrexTest::SetUp();
        {
            amrex::ParmParse pp("Actuator.FixedWingLine");
            const amrex::Vector<amrex::Real> eps{2.0_rt, 2.0_rt, 2.0_rt};
            const amrex::Vector<amrex::Real> eps_chord{
                0.25_rt, 0.25_rt, 0.25_rt};
            pp.addarr("epsilon", eps);
            pp.addarr("epsilon_chord", eps_chord);
        }
    }
};

//! Message of the abort of fllc_init, empty when it does not abort
std::string fllc_init_abort_message(
    FLLCData& data, const ComponentView& view, const amrex::Real eps_chord)
{
    try {
        fllc_init(data, view, eps_chord);
    } catch (const amrex::RuntimeError& err) {
        return err.what();
    }
    return "";
}

} // namespace

TEST_F(FLLCParseTest, reads_inputs_in_range)
{
    amrex::ParmParse pp("Actuator.FixedWingLine");
    pp.add("fllc_relaxation_factor", 1.0_rt);
    pp.add("fllc_epsilon_dr_ratio", 3.0_rt);
    const utils::ActParser ap("Actuator.FixedWingLine", "Actuator.F1");
    FLLCData data;

    fllc_parse(ap, data);

    EXPECT_EQ(data.epsilon, 2.0_rt);
    EXPECT_EQ(data.relaxation_factor, 1.0_rt);
    EXPECT_EQ(data.eps_dr, 3.0_rt);
}

TEST_F(FLLCParseTest, rejects_relaxation_factor_of_zero)
{
    // With a factor of 0 the correction stays at zero, without a word
    amrex::ParmParse pp("Actuator.FixedWingLine");
    pp.add("fllc_relaxation_factor", 0.0_rt);
    const utils::ActParser ap("Actuator.FixedWingLine", "Actuator.F1");
    FLLCData data;

    EXPECT_THROW(fllc_parse(ap, data), amrex::RuntimeError);
}

TEST_F(FLLCParseTest, rejects_relaxation_factor_above_one)
{
    // Above 1 every update of the correction overshoots
    amrex::ParmParse pp("Actuator.FixedWingLine");
    pp.add("fllc_relaxation_factor", 2.0_rt);
    const utils::ActParser ap("Actuator.FixedWingLine", "Actuator.F1");
    FLLCData data;

    EXPECT_THROW(fllc_parse(ap, data), amrex::RuntimeError);
}

TEST_F(FLLCParseTest, rejects_epsilon_dr_ratio_of_zero)
{
    // The non-uniform points are spaced by the optimal epsilon over this
    // ratio
    amrex::ParmParse pp("Actuator.FixedWingLine");
    pp.add("fllc_epsilon_dr_ratio", 0.0_rt);
    const utils::ActParser ap("Actuator.FixedWingLine", "Actuator.F1");
    FLLCData data;

    EXPECT_THROW(fllc_parse(ap, data), amrex::RuntimeError);
}

TEST(TestFLLCData, init_rejects_zero_chord)
{
    const int num_points = 5;

    // The kernel of the correction divides by chord times epsilon_chord
    FLLCTestBlade blade(num_points, num_points, vs::Vector::zero());
    blade.chord[2] = 0.0_rt;
    FLLCData data;
    data.epsilon = 2.0_rt;
    data.nonuniform = false;

    EXPECT_THROW(fllc_init(data, blade.view, 1.0_rt), amrex::RuntimeError);
}

TEST(TestFLLCData, init_names_zero_chord_with_nonuniform_points)
{
    const int num_points = 5;

    FLLCTestBlade blade(num_points, num_points, vs::Vector::zero());
    blade.chord[2] = 0.0_rt;
    FLLCData data;
    data.epsilon = 2.0_rt;
    data.nonuniform = true;
    data.eps_dr = 1.0_rt;

    // The message has to tell what is wrong with the inputs
    const auto msg = fllc_init_abort_message(data, blade.view, 1.0_rt);
    EXPECT_NE(msg.find("chord must be larger than 0"), std::string::npos)
        << "message: " << msg;
}

TEST(TestFLLCData, init_rejects_zero_epsilon_chord)
{
    const int num_points = 5;

    FLLCTestBlade blade(num_points, num_points, vs::Vector::zero());
    FLLCData data;
    data.epsilon = 2.0_rt;
    data.nonuniform = false;

    EXPECT_THROW(fllc_init(data, blade.view, 0.0_rt), amrex::RuntimeError);
}

TEST(TestFLLCData, init_rejects_coincident_points)
{
    const int num_points = 5;

    // Two actuator points at the same place: the line has no length there
    FLLCTestBlade blade(num_points, num_points, vs::Vector::zero());
    blade.pos[3] = blade.pos[2];
    FLLCData data;
    data.epsilon = 2.0_rt;
    data.nonuniform = false;

    EXPECT_THROW(fllc_init(data, blade.view, 1.0_rt), amrex::RuntimeError);
}

TEST(TestFLLCData, init_stops_when_nonuniform_points_cannot_reach_tip)
{
    const int num_points = 5;

    // A tip chord of almost zero: the non-uniform points get a spacing
    // there that never brings them to the tip
    FLLCTestBlade blade(num_points, num_points, vs::Vector::zero());
    blade.chord[num_points - 1] = 3.0e-16_rt;
    FLLCData data;
    data.epsilon = 2.0_rt;
    data.nonuniform = true;
    data.eps_dr = 1.0_rt;

    testing::internal::CaptureStdout();
    const auto msg = fllc_init_abort_message(data, blade.view, 1.0_rt);
    testing::internal::GetCapturedStdout();
    EXPECT_NE(
        msg.find("cannot place the non-uniform points"), std::string::npos)
        << "message: " << msg;
}

TEST(TestFLLCData, init_warns_when_points_do_not_resolve_optimal_epsilon)
{
    const int num_points = 5;

    // Actuator points one unit apart, and an optimal epsilon of a quarter
    // of that
    FLLCTestBlade blade(num_points, num_points, vs::Vector::zero());
    FLLCData data;
    data.epsilon = 2.0_rt;
    data.nonuniform = true;
    data.eps_dr = 1.0_rt;

    testing::internal::CaptureStdout();
    fllc_init(data, blade.view, 0.25_rt);
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_NE(out.find("WARNING: FLLC"), std::string::npos);
    EXPECT_NE(out.find("at 5 of 5 points"), std::string::npos);
}

TEST(TestFLLCData, init_is_silent_when_points_resolve_optimal_epsilon)
{
    const int num_points = 5;

    // Actuator points one unit apart, and an optimal epsilon of twice that
    FLLCTestBlade blade(num_points, num_points, vs::Vector::zero());
    FLLCData data;
    data.epsilon = 2.0_rt;
    data.nonuniform = true;
    data.eps_dr = 1.0_rt;

    testing::internal::CaptureStdout();
    fllc_init(data, blade.view, 2.0_rt);
    const std::string out = testing::internal::GetCapturedStdout();
    EXPECT_EQ(out.find("WARNING"), std::string::npos);
}

} // namespace kynema_sgf::actuator
