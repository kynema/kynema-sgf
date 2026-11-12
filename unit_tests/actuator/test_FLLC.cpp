#include "ks_test_utils/AmrexTest.H"
#include "src/actuator/FLLC.H"
#include "src/actuator/FLLCOp.H"
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

} // namespace kynema_sgf::actuator
