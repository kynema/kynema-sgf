#include "aw_test_utils/MeshTest.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/fvm/laplacian.H"
#include "amr-wind/fvm/divergence.H"
#include "AnalyticalFunction.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

class FvmOperatorsTest : public MeshTest
{};

namespace {

void initialize_velocity(
    const amrex::Geometry& geom,
    const amrex::Box& bx,
    const int pdegree,
    amrex::Gpu::DeviceVector<amrex::Real>& cu,
    amrex::Gpu::DeviceVector<amrex::Real>& cv,
    amrex::Gpu::DeviceVector<amrex::Real>& cw,
    const amrex::Array4<amrex::Real>& vel_arr)
{
    auto problo = geom.ProbLoArray();
    auto probhi = geom.ProbHiArray();
    auto dx = geom.CellSizeArray();

    const amrex::Real* cu_ptr = cu.data();
    const amrex::Real* cv_ptr = cv.data();
    const amrex::Real* cw_ptr = cw.data();

    // grow the box by 1 so that x,y,z go out of bounds and min(max()) corrects
    // it and it fills the ghosts with wall values
    amrex::ParallelFor(grow(bx, 1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real x = amrex::min(
            amrex::max<amrex::Real>(
                problo[0] + ((i + 0.5_rt) * dx[0]), problo[0]),
            probhi[0]);
        const amrex::Real y = amrex::min(
            amrex::max<amrex::Real>(
                problo[1] + ((j + 0.5_rt) * dx[1]), problo[1]),
            probhi[1]);
        const amrex::Real z = amrex::min(
            amrex::max<amrex::Real>(
                problo[2] + ((k + 0.5_rt) * dx[2]), problo[2]),
            probhi[2]);

        vel_arr(i, j, k, 0) =
            analytical_function::phi_eval(pdegree, cu_ptr, x, y, z);
        vel_arr(i, j, k, 1) =
            analytical_function::phi_eval(pdegree, cv_ptr, x, y, z);
        vel_arr(i, j, k, 2) =
            analytical_function::phi_eval(pdegree, cw_ptr, x, y, z);
    });
}

amrex::Real grad_test_impl(amr_wind::Field& vel, const int pdegree)
{

    const int ncoeff = (pdegree + 1) * (pdegree + 1) * (pdegree + 1);

    amrex::Gpu::DeviceVector<amrex::Real> cu(ncoeff, 0.00123_rt);
    amrex::Gpu::DeviceVector<amrex::Real> cv(ncoeff, 0.00213_rt);
    amrex::Gpu::DeviceVector<amrex::Real> cw(ncoeff, 0.00346_rt);

    const auto& geom = vel.repo().mesh().Geom();

    run_algorithm(vel, [&](const int lev, const amrex::MFIter& mfi) {
        auto vel_arr = vel(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_velocity(geom[lev], bx, pdegree, cu, cv, cw, vel_arr);
    });

    auto grad_vel = amr_wind::fvm::gradient(vel);

    const int nlevels = vel.repo().num_active_levels();
    amrex::Real error_total = 0.0_rt;

    const amrex::Real* cu_ptr = cu.data();
    const amrex::Real* cv_ptr = cv.data();
    const amrex::Real* cw_ptr = cw.data();

    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& problo = geom[lev].ProbLoArray();
        const auto& dx = geom[lev].CellSizeArray();

        error_total += amrex::ReduceSum(
            (*grad_vel)(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& gvel) -> amrex::Real {
                amrex::Real error = 0.0_rt;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                    const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
                    const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

                    error += std::abs(
                        gvel(i, j, k, 0) - analytical_function::dphidx_eval(
                                               pdegree, cu_ptr, x, y, z));
                    error += std::abs(
                        gvel(i, j, k, 1) - analytical_function::dphidy_eval(
                                               pdegree, cu_ptr, x, y, z));
                    error += std::abs(
                        gvel(i, j, k, 2) - analytical_function::dphidz_eval(
                                               pdegree, cu_ptr, x, y, z));

                    error += std::abs(
                        gvel(i, j, k, 3) - analytical_function::dphidx_eval(
                                               pdegree, cv_ptr, x, y, z));
                    error += std::abs(
                        gvel(i, j, k, 4) - analytical_function::dphidy_eval(
                                               pdegree, cv_ptr, x, y, z));
                    error += std::abs(
                        gvel(i, j, k, 5) - analytical_function::dphidz_eval(
                                               pdegree, cv_ptr, x, y, z));

                    error += std::abs(
                        gvel(i, j, k, 6) - analytical_function::dphidx_eval(
                                               pdegree, cw_ptr, x, y, z));
                    error += std::abs(
                        gvel(i, j, k, 7) - analytical_function::dphidy_eval(
                                               pdegree, cw_ptr, x, y, z));
                    error += std::abs(
                        gvel(i, j, k, 8) - analytical_function::dphidz_eval(
                                               pdegree, cw_ptr, x, y, z));
                });

                return error;
            });
    }

    return error_total;
}

amrex::Real laplacian_test_impl(amr_wind::Field& vel, const int pdegree)
{

    const int ncoeff = (pdegree + 1) * (pdegree + 1) * (pdegree + 1);

    amrex::Gpu::DeviceVector<amrex::Real> cu(ncoeff, 0.00123_rt);
    amrex::Gpu::DeviceVector<amrex::Real> cv(ncoeff, 0.00213_rt);
    amrex::Gpu::DeviceVector<amrex::Real> cw(ncoeff, 0.00346_rt);

    const auto& geom = vel.repo().mesh().Geom();

    run_algorithm(vel, [&](const int lev, const amrex::MFIter& mfi) {
        auto vel_arr = vel(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_velocity(geom[lev], bx, pdegree, cu, cv, cw, vel_arr);
    });

    auto lap = amr_wind::fvm::laplacian(vel);

    const int nlevels = vel.repo().num_active_levels();
    amrex::Real error_total = 0.0_rt;

    const amrex::Real* cu_ptr = cu.data();
    const amrex::Real* cv_ptr = cv.data();
    const amrex::Real* cw_ptr = cw.data();

    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& problo = geom[lev].ProbLoArray();
        const auto& dx = geom[lev].CellSizeArray();

        error_total += amrex::ReduceSum(
            (*lap)(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& lap_arr)
                -> amrex::Real {
                amrex::Real error = 0.0_rt;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                    const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
                    const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

                    error += std::abs(
                        lap_arr(i, j, k) -
                        analytical_function::laplacian(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                });

                return error;
            });
    }

    return error_total;
}

amrex::Real divergence_test_impl(amr_wind::Field& vel, const int pdegree)
{

    const int ncoeff = (pdegree + 1) * (pdegree + 1) * (pdegree + 1);

    amrex::Gpu::DeviceVector<amrex::Real> cu(ncoeff, 0.00123_rt);
    amrex::Gpu::DeviceVector<amrex::Real> cv(ncoeff, 0.00213_rt);
    amrex::Gpu::DeviceVector<amrex::Real> cw(ncoeff, 0.00346_rt);

    const auto& geom = vel.repo().mesh().Geom();

    run_algorithm(vel, [&](const int lev, const amrex::MFIter& mfi) {
        auto vel_arr = vel(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_velocity(geom[lev], bx, pdegree, cu, cv, cw, vel_arr);
    });

    auto div = amr_wind::fvm::divergence(vel);

    const int nlevels = vel.repo().num_active_levels();
    amrex::Real error_total = 0.0_rt;

    const amrex::Real* cu_ptr = cu.data();
    const amrex::Real* cv_ptr = cv.data();
    const amrex::Real* cw_ptr = cw.data();

    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& problo = geom[lev].ProbLoArray();
        const auto& dx = geom[lev].CellSizeArray();

        error_total += amrex::ReduceSum(
            (*div)(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& div_arr)
                -> amrex::Real {
                amrex::Real error = 0.0_rt;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                    const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
                    const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

                    error += std::abs(
                        div_arr(i, j, k) -
                        analytical_function::divergence(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                });

                return error;
            });
    }

    return error_total;
}

} // namespace

TEST_F(FvmOperatorsTest, gradient)
{

    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt;

    populate_parameters();
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<int> periodic{{0, 0, 0}};
        pp.addarr("is_periodic", periodic);
    }

    initialize_mesh();

    auto& repo = sim().repo();
    const int ncomp = 3;
    const int nghost = 1;
    auto& vel = repo.declare_field("vel", ncomp, nghost);
    const int pdegree = 2;
    auto error_total = grad_test_impl(vel, pdegree);

    amrex::ParallelDescriptor::ReduceRealSum(error_total);

    EXPECT_NEAR(error_total, 0.0_rt, tol);
}

TEST_F(FvmOperatorsTest, laplacian)
{

    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e5_rt;

    populate_parameters();
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<int> periodic{{0, 0, 0}};
        pp.addarr("is_periodic", periodic);
    }

    initialize_mesh();

    auto& repo = sim().repo();
    const int ncomp = 3;
    const int nghost = 1;
    auto& vel = repo.declare_field("vel", ncomp, nghost);

    const int pdegree = 2;
    auto error_total = laplacian_test_impl(vel, pdegree);

    amrex::ParallelDescriptor::ReduceRealSum(error_total);

    EXPECT_NEAR(error_total, 0.0_rt, tol);
}

TEST_F(FvmOperatorsTest, divergence)
{

    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e5_rt;

    populate_parameters();
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<int> periodic{{0, 0, 0}};
        pp.addarr("is_periodic", periodic);
    }

    initialize_mesh();

    auto& repo = sim().repo();
    const int ncomp = 3;
    const int nghost = 1;
    auto& vel = repo.declare_field("vel", ncomp, nghost);

    const int pdegree = 2;
    auto error_total = divergence_test_impl(vel, pdegree);

    amrex::ParallelDescriptor::ReduceRealSum(error_total);

    EXPECT_NEAR(error_total, 0.0_rt, tol);
}

} // namespace amr_wind_tests
