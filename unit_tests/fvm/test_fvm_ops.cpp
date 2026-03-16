#include "aw_test_utils/MeshTest.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/fvm/strainrate.H"
#include "amr-wind/fvm/vorticity.H"
#include "amr-wind/fvm/vorticity_mag.H"
#include "amr-wind/fvm/qcriterion.H"
#include "amr-wind/fvm/laplacian.H"
#include "amr-wind/fvm/divergence.H"
#include "amr-wind/fvm/curvature.H"
#include "amr-wind/fvm/nonLinearSum.H"
#include "AnalyticalFunction.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

class FvmOpTest : public MeshTest
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

amrex::Real strainrate_test_impl(amr_wind::Field& vel, const int pdegree)
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

    auto str = amr_wind::fvm::strainrate(vel);

    const int nlevels = vel.repo().num_active_levels();
    amrex::Real error_total = 0.0_rt;

    const amrex::Real* cu_ptr = cu.data();
    const amrex::Real* cv_ptr = cv.data();
    const amrex::Real* cw_ptr = cw.data();

    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& problo = geom[lev].ProbLoArray();
        const auto& dx = geom[lev].CellSizeArray();

        error_total += amrex::ReduceSum(
            (*str)(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& str_arr)
                -> amrex::Real {
                amrex::Real error = 0.0_rt;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                    const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
                    const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

                    error += std::abs(
                        str_arr(i, j, k) -
                        analytical_function::strainrate(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                });

                return error;
            });
    }

    return error_total;
}

amrex::Real nonlinearsum_test_impl(amr_wind::Field& vel, const int pdegree)
{

    const int ncoeff = (pdegree + 1) * (pdegree + 1) * (pdegree + 1);

    amrex::Gpu::DeviceVector<amrex::Real> cu(ncoeff, 0.000193_rt);
    amrex::Gpu::DeviceVector<amrex::Real> cv(ncoeff, 0.000463_rt);
    amrex::Gpu::DeviceVector<amrex::Real> cw(ncoeff, 0.000386_rt);

    const auto& geom = vel.repo().mesh().Geom();

    run_algorithm(vel, [&](const int lev, const amrex::MFIter& mfi) {
        auto vel_arr = vel(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_velocity(geom[lev], bx, pdegree, cu, cv, cw, vel_arr);
    });

    auto nonlinearsum = amr_wind::fvm::nonlinearsum(vel);

    const int nlevels = vel.repo().num_active_levels();
    amrex::Real error_total = 0.0_rt;

    const amrex::Real* cu_ptr = cu.data();
    const amrex::Real* cv_ptr = cv.data();
    const amrex::Real* cw_ptr = cw.data();

    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& problo = geom[lev].ProbLoArray();
        const auto& dx = geom[lev].CellSizeArray();

        error_total += amrex::ReduceSum(
            (*nonlinearsum)(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& nonlinearsum_arr)
                -> amrex::Real {
                amrex::Real error = 0.0_rt;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                    const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
                    const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

                    error += std::abs(
                        nonlinearsum_arr(i, j, k, 0) -
                        analytical_function::nonlinearsum_11(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                    error += std::abs(
                        nonlinearsum_arr(i, j, k, 1) -
                        analytical_function::nonlinearsum_12(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                    error += std::abs(
                        nonlinearsum_arr(i, j, k, 2) -
                        analytical_function::nonlinearsum_13(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                    error += std::abs(
                        nonlinearsum_arr(i, j, k, 3) -
                        analytical_function::nonlinearsum_21(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                    error += std::abs(
                        nonlinearsum_arr(i, j, k, 4) -
                        analytical_function::nonlinearsum_22(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                    error += std::abs(
                        nonlinearsum_arr(i, j, k, 5) -
                        analytical_function::nonlinearsum_23(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                    error += std::abs(
                        nonlinearsum_arr(i, j, k, 6) -
                        analytical_function::nonlinearsum_31(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                    error += std::abs(
                        nonlinearsum_arr(i, j, k, 7) -
                        analytical_function::nonlinearsum_32(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                    error += std::abs(
                        nonlinearsum_arr(i, j, k, 8) -
                        analytical_function::nonlinearsum_33(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                });

                return error;
            });
    }
    return error_total;
}

amrex::Real vorticity_test_impl(amr_wind::Field& vel, const int pdegree)
{

    const int ncoeff = (pdegree + 1) * (pdegree + 1) * (pdegree + 1);

    amrex::Gpu::DeviceVector<amrex::Real> cu(ncoeff, 0.000193_rt);
    amrex::Gpu::DeviceVector<amrex::Real> cv(ncoeff, 0.000463_rt);
    amrex::Gpu::DeviceVector<amrex::Real> cw(ncoeff, 0.000386_rt);

    const auto& geom = vel.repo().mesh().Geom();

    run_algorithm(vel, [&](const int lev, const amrex::MFIter& mfi) {
        auto vel_arr = vel(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_velocity(geom[lev], bx, pdegree, cu, cv, cw, vel_arr);
    });

    auto vorticity = amr_wind::fvm::vorticity(vel);

    const int nlevels = vel.repo().num_active_levels();
    amrex::Real error_total = 0.0_rt;

    const amrex::Real* cu_ptr = cu.data();
    const amrex::Real* cv_ptr = cv.data();
    const amrex::Real* cw_ptr = cw.data();

    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& problo = geom[lev].ProbLoArray();
        const auto& dx = geom[lev].CellSizeArray();

        error_total += amrex::ReduceSum(
            (*vorticity)(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& vort_arr)
                -> amrex::Real {
                amrex::Real error = 0.0_rt;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                    const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
                    const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

                    error += std::abs(
                        vort_arr(i, j, k, 0) -
                        analytical_function::vorticity_x(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                    error += std::abs(
                        vort_arr(i, j, k, 1) -
                        analytical_function::vorticity_y(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                    error += std::abs(
                        vort_arr(i, j, k, 2) -
                        analytical_function::vorticity_z(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                });

                return error;
            });
    }
    return error_total;
}

amrex::Real vorticity_mag_test_impl(amr_wind::Field& vel, const int pdegree)
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

    auto vrt_mag = amr_wind::fvm::vorticity_mag(vel);
    auto vorticity = amr_wind::fvm::vorticity(vel);

    const int nlevels = vel.repo().num_active_levels();
    amrex::Real error_total = 0.0_rt;

    const amrex::Real* cu_ptr = cu.data();
    const amrex::Real* cv_ptr = cv.data();
    const amrex::Real* cw_ptr = cw.data();

    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& problo = geom[lev].ProbLoArray();
        const auto& dx = geom[lev].CellSizeArray();

        error_total += amrex::ReduceSum(
            (*vorticity)(lev), (*vrt_mag)(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& vrt_arr,
                amrex::Array4<amrex::Real const> const& vrt_mag_arr)
                -> amrex::Real {
                amrex::Real error = 0.0_rt;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                    const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
                    const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

                    error += std::abs(
                        vrt_mag_arr(i, j, k) -
                        analytical_function::vorticity_mag(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));

                    const amrex::Real vortmag = std::sqrt(
                        (vrt_arr(i, j, k, 0) * vrt_arr(i, j, k, 0)) +
                        (vrt_arr(i, j, k, 1) * vrt_arr(i, j, k, 1)) +
                        (vrt_arr(i, j, k, 2) * vrt_arr(i, j, k, 2)));

                    error += std::abs(vrt_mag_arr(i, j, k) - vortmag);
                });

                return error;
            });
    }

    return error_total;
}

amrex::Real q_criterion_test_impl(amr_wind::Field& vel, const int pdegree)
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

    auto qcrit = amr_wind::fvm::q_criterion(vel);

    const int nlevels = vel.repo().num_active_levels();
    amrex::Real error_total = 0.0_rt;

    const amrex::Real* cu_ptr = cu.data();
    const amrex::Real* cv_ptr = cv.data();
    const amrex::Real* cw_ptr = cw.data();

    for (int lev = 0; lev < nlevels; ++lev) {

        const auto& problo = geom[lev].ProbLoArray();
        const auto& dx = geom[lev].CellSizeArray();

        error_total += amrex::ReduceSum(
            (*qcrit)(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& qcrit_arr)
                -> amrex::Real {
                amrex::Real error = 0.0_rt;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                    const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
                    const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);

                    error += std::abs(
                        qcrit_arr(i, j, k) -
                        analytical_function::q_criterion(
                            pdegree, cu_ptr, cv_ptr, cw_ptr, x, y, z));
                });

                return error;
            });
    }

    return error_total;
}

} // namespace

TEST_F(FvmOpTest, nonlinearsum)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e7_rt;

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
    auto error_total = nonlinearsum_test_impl(vel, pdegree);

    amrex::ParallelDescriptor::ReduceRealSum(error_total);

    EXPECT_NEAR(error_total, 0.0_rt, tol);
}

TEST_F(FvmOpTest, strainrate)
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
    auto error_total = strainrate_test_impl(vel, pdegree);

    amrex::ParallelDescriptor::ReduceRealSum(error_total);

    EXPECT_NEAR(error_total, 0.0_rt, tol);
}

TEST_F(FvmOpTest, vorticity)
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
    auto error_total = vorticity_test_impl(vel, pdegree);

    amrex::ParallelDescriptor::ReduceRealSum(error_total);

    EXPECT_NEAR(error_total, 0.0_rt, tol);
}

TEST_F(FvmOpTest, vorticity_mag)
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
    auto error_total = vorticity_mag_test_impl(vel, pdegree);

    amrex::ParallelDescriptor::ReduceRealSum(error_total);

    EXPECT_NEAR(error_total, 0.0_rt, tol);
}

TEST_F(FvmOpTest, q_criterion)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e7_rt;

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
    auto error_total = q_criterion_test_impl(vel, pdegree);

    amrex::ParallelDescriptor::ReduceRealSum(error_total);

    EXPECT_NEAR(error_total, 0.0_rt, tol);
}

} // namespace amr_wind_tests
