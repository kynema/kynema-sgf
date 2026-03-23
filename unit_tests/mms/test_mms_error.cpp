#include "mms_test_utils.H"
#include "amr-wind/utilities/trig_ops.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"

#include "masa.h"
#include "amr-wind/physics/mms/MMS.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

namespace {
void perturb_vel_field(
    const amrex::Box& bx,
    const amrex::Array4<amrex::Real>& vel,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& problo)
{
    amrex::ParallelFor(
        bx, [vel, dx, problo] AMREX_GPU_DEVICE(int i, int j, int k) {
            const amrex::Real x = problo[0] + (i + 0.5_rt) * dx[0];
            const amrex::Real y = problo[1] + (j + 0.5_rt) * dx[1];
            const amrex::Real z = problo[2] + (k + 0.5_rt) * dx[2];
            vel(i, j, k, 0) += std::sin(amr_wind::utils::two_pi() * x);
            vel(i, j, k, 1) += std::sin(amr_wind::utils::two_pi() * 2.0_rt * y);
            vel(i, j, k, 2) += std::cos(amr_wind::utils::two_pi() * z);
        });
}
} // namespace

TEST_F(MMSMeshTest, mms_error)
{
#if defined(AMREX_USE_HIP)
    GTEST_SKIP();
#else
    populate_parameters();
    utils::populate_mms_params();

    initialize_mesh();
    auto& frepo = mesh().field_repo();
    auto& velocityf = frepo.declare_field("velocity", 3, 0);
    frepo.declare_field("density");

    amr_wind::mms::MMS mms(sim());
    const int nlevels = mesh().num_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        mms.initialize_fields(lev, mesh().Geom(lev));
    }

    // perturb the velocity fields from the initial condition
    auto velocity = velocityf.vec_ptrs();
    run_algorithm(
        mesh().num_levels(), velocity,
        [&](const int lev, const amrex::MFIter& mfi) {
            auto vel = velocity[lev]->array(mfi);

            const auto& bx = mfi.validbox();
            const auto& dx = mesh().Geom(lev).CellSizeArray();
            const auto& problo = mesh().Geom(lev).ProbLoArray();

            perturb_vel_field(bx, vel, dx, problo);
        });

    const amrex::Real u_mms_err =
        mms.compute_error(0, velocityf, masa_eval_3d_exact_u);
    const amrex::Real v_mms_err =
        mms.compute_error(1, velocityf, masa_eval_3d_exact_v);
    const amrex::Real w_mms_err =
        mms.compute_error(2, velocityf, masa_eval_3d_exact_w);

    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    amrex::Array<amrex::Real, AMREX_SPACEDIM> golds = {
        0.67158428586284569_rt, 0.6978702996158761_rt, 0.74092816587175314_rt};
    EXPECT_NEAR(u_mms_err, golds[0], tol);
    EXPECT_NEAR(v_mms_err, golds[1], tol);
    EXPECT_NEAR(w_mms_err, golds[2], tol);
#endif
}

} // namespace amr_wind_tests
