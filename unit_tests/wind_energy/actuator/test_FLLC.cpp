#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/wind_energy/actuator/FLLC.H"
#include "amr-wind/wind_energy/actuator/FLLCOp.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator {

TEST(TestFLLCData, data_initializes_with_cviews)
{
    const int num_points = 10;

    VecList dummy_vec(num_points);
    RealList dummy_real(num_points);
    TensorList dummy_tensor(num_points);

    FLLCData data;

    ComponentView view;
    view.pos = ::amr_wind::utils::slice(dummy_vec, 0, num_points);
    view.force = ::amr_wind::utils::slice(dummy_vec, 0, num_points);
    view.epsilon = ::amr_wind::utils::slice(dummy_vec, 0, num_points);
    view.orientation = ::amr_wind::utils::slice(dummy_tensor, 0, num_points);
    view.chord = ::amr_wind::utils::slice(dummy_real, 0, num_points);
    view.vel_rel = ::amr_wind::utils::slice(dummy_vec, 0, num_points);

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
    view.pos = ::amr_wind::utils::slice(dummy_vec, 0, num_points);
    view.force = ::amr_wind::utils::slice(dummy_vec, 0, num_points);
    view.epsilon = ::amr_wind::utils::slice(dummy_vec, 0, num_points);
    view.orientation = ::amr_wind::utils::slice(dummy_tensor, 0, num_points);
    view.chord = ::amr_wind::utils::slice(dummy_real, 0, num_points);
    view.vel_rel = ::amr_wind::utils::slice(dummy_vec, 0, num_points);

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

} // namespace amr_wind::actuator
