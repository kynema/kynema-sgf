#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

class VOFOpTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{4, 4, 4}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 4);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{{1.0_rt, 1.0_rt, 1.0_rt}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

namespace {

void initialize_volume_fractions(
    const int dir,
    const amrex::Box& bx,
    const amrex::Array4<amrex::Real>& vof_arr)
{
    // grow the box by 1 so that x,y,z go out of bounds and min(max()) corrects
    // it and it fills the ghosts with wall values
    amrex::ParallelFor(grow(bx, 1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        int ii = (dir != 0 ? i : 0);
        int jj = (dir != 1 ? j : 0);
        int kk = (dir != 2 ? k : 0);
        if (ii + jj + kk > 3) {
            vof_arr(i, j, k) = 0.0_rt;
        }
        if (ii + jj + kk == 3) {
            vof_arr(i, j, k) = 0.5_rt;
        }
        if (ii + jj + kk < 3) {
            vof_arr(i, j, k) = 1.0_rt;
        }
    });
}

void init_vof(amr_wind::Field& vof, const int dir)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto vof_arr = vof(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_volume_fractions(dir, bx, vof_arr);
    });
}

void initialize_volume_fractions_horizontal(
    const int dir,
    const amrex::Box& bx,
    const amrex::Real vof_val,
    const amrex::Array4<amrex::Real>& vof_arr)
{
    // grow the box by 1 so that x,y,z go out of bounds and min(max()) corrects
    // it and it fills the ghosts with wall values
    const amrex::Real vv = vof_val;
    amrex::ParallelFor(grow(bx, 1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        int ii = (dir == 0 ? i : (dir == 1 ? j : k));
        if (ii > 1) {
            vof_arr(i, j, k) = 0.0_rt;
        }
        if (ii == 1) {
            vof_arr(i, j, k) = vv;
        }
        if (ii < 1) {
            vof_arr(i, j, k) = 1.0_rt;
        }
    });
}

void init_vof_h(amr_wind::Field& vof, const amrex::Real vof_val, const int dir)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto vof_arr = vof(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_volume_fractions_horizontal(dir, bx, vof_val, vof_arr);
    });
}

void initialize_volume_fractions(
    const amrex::Box& bx, const amrex::Array4<amrex::Real>& vof_arr)
{
    amrex::ParallelFor(grow(bx, 1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        vof_arr(i, j, k) =
            (0.13_rt * (static_cast<amrex::Real>(i) - 1.5_rt)) +
            (0.04_rt * std::pow(static_cast<amrex::Real>(j) - 1.0_rt, 2.0_rt)) +
            (0.01_rt * std::pow(static_cast<amrex::Real>(k) - 2.0_rt, 3.0_rt)) +
            0.5_rt;
    });
}

void init_vof(amr_wind::Field& vof)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto vof_arr = vof(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_volume_fractions(bx, vof_arr);
    });
}

void initialize_iblank_distribution(
    const int dir, const amrex::Box& bx, const amrex::Array4<int>& iblk_arr)
{
    amrex::ParallelFor(grow(bx, 1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        int ii = (dir != 0 ? i : 0);
        int jj = (dir != 1 ? j : 0);
        int kk = (dir != 2 ? k : 0);
        if (ii + jj + kk > 4 || ii + jj + kk < 2) {
            iblk_arr(i, j, k) = -1;
        } else {
            iblk_arr(i, j, k) = 1;
        }
    });
}

void init_iblank(amr_wind::IntField& iblank, const int dir)
{
    run_algorithm(iblank, [&](const int lev, const amrex::MFIter& mfi) {
        auto iblk_arr = iblank(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_iblank_distribution(dir, bx, iblk_arr);
    });
}

void initialize_volume_fractions_iblank(
    const int dir,
    const amrex::Box& bx,
    const amrex::Array4<amrex::Real>& vof_arr,
    const amrex::Array4<int>& iblk_arr)
{
    // Does a horizontal interface that is different outside of iblank region
    amrex::ParallelFor(grow(bx, 1), [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // Turns on index tangent to interface
        int iit = (dir != 0 ? i : 1);
        int jjt = (dir != 1 ? j : 1);
        int kkt = (dir != 2 ? k : 1);
        // Turns on index normal to interface
        int iin = (dir == 0 ? i : 0);
        int jjn = (dir == 1 ? j : 0);
        int kkn = (dir == 2 ? k : 0);
        // Ordinary vof distribution
        if (iin + jjn + kkn > 2) {
            vof_arr(i, j, k) = 0.0_rt;
        }
        if (iin + jjn + kkn == 2) {
            vof_arr(i, j, k) = 0.5_rt;
        }
        if (iin + jjn + kkn < 2) {
            vof_arr(i, j, k) = 1.0_rt;
        }
        // iblank distribution
        iblk_arr(i, j, k) = 1;
        if (iit > 2 || iit < 1 || jjt > 2 || jjt < 1 || kkt > 2 || kkt < 1) {
            // Restrict iblank = 1 to central block
            iblk_arr(i, j, k) = -1;
            // Change vof where iblank = -1
            if (iin + jjn + kkn == 3) {
                vof_arr(i, j, k) = 0.25_rt;
            }
            if (iin + jjn + kkn == 1) {
                vof_arr(i, j, k) = 0.75_rt;
            }
        }
    });
}

void init_vof_iblank(
    amr_wind::Field& vof, amr_wind::IntField& iblank, const int dir)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto vof_arr = vof(lev).array(mfi);
        auto iblk_arr = iblank(lev).array(mfi);
        const auto& bx = mfi.validbox();
        initialize_volume_fractions_iblank(dir, bx, vof_arr, iblk_arr);
    });
}

amrex::Real normal_vector_test_impl(amr_wind::Field& vof, const int dir)
{
    amrex::Real error_total = 0.0_rt;

    for (int lev = 0; lev < vof.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            vof(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& vof_arr)
                -> amrex::Real {
                amrex::Real error = 0.0_rt;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    amrex::Real mx, my, mz;
                    amr_wind::multiphase::mixed_youngs_central_normal(
                        i, j, k, vof_arr, mx, my, mz);

                    int ii = (dir != 0 ? i : 0);
                    int jj = (dir != 1 ? j : 0);
                    int kk = (dir != 2 ? k : 0);

                    // Use L1 norm, check cells where slope is known
                    if (ii + jj + kk == 3) {
                        error += std::abs(mx - (dir != 0 ? 0.5_rt : 0.0_rt));
                        error += std::abs(my - (dir != 1 ? 0.5_rt : 0.0_rt));
                        error += std::abs(mz - (dir != 2 ? 0.5_rt : 0.0_rt));
                    }
                });

                return error;
            });
    }
    return error_total;
}

amrex::Real normal_vector_neumann_test_impl(
    amr_wind::Field& vof, amr_wind::IntField& iblk_fld)
{
    amrex::Real error_total = 0.0_rt;

    for (int lev = 0; lev < vof.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            vof(lev), iblk_fld(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& vof_arr,
                amrex::Array4<int const> const& iblank) -> amrex::Real {
                amrex::Real error = 0.0_rt;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    int ibdy =
                        (iblank(i, j, k) != iblank(i - 1, j, k)) ? -1 : 0;
                    int jbdy =
                        (iblank(i, j, k) != iblank(i, j - 1, k)) ? -1 : 0;
                    int kbdy =
                        (iblank(i, j, k) != iblank(i, j, k - 1)) ? -1 : 0;
                    ibdy = (iblank(i, j, k) != iblank(i + 1, j, k)) ? +1 : ibdy;
                    jbdy = (iblank(i, j, k) != iblank(i, j + 1, k)) ? +1 : jbdy;
                    kbdy = (iblank(i, j, k) != iblank(i, j, k + 1)) ? +1 : kbdy;
                    amrex::Real mxn, myn, mzn;
                    amr_wind::multiphase::
                        youngs_finite_difference_normal_neumann(
                            i, j, k, ibdy, jbdy, kbdy, vof_arr, mxn, myn, mzn);
                    amrex::Real mx, my, mz;
                    amr_wind::multiphase::youngs_finite_difference_normal(
                        i, j, k, vof_arr, mx, my, mz);

                    // Use L1 norm, check against non-neumann implementation
                    // Slope across overset boundary should be different
                    constexpr amrex::Real slp_tol =
                        std::numeric_limits<float>::epsilon() * 1.0e1_rt;
                    if (ibdy != 0) {
                        // x slope should be different
                        error += std::abs(mx - mxn) > slp_tol ? 0.0_rt : 1.0_rt;
                    }
                    if (jbdy != 0) {
                        // y slope should be different
                        error += std::abs(my - myn) > slp_tol ? 0.0_rt : 1.0_rt;
                    }
                    if (kbdy != 0) {
                        // z slope should be different
                        error += std::abs(mz - mzn) > slp_tol ? 0.0_rt : 1.0_rt;
                    }
                    // Slope should otherwise be the same
                    if (ibdy == 0 && jbdy == 0 && kbdy == 0) {
                        error += std::abs(mx - mxn);
                        error += std::abs(my - myn);
                        error += std::abs(mz - mzn);
                    }
                });

                return error;
            });
    }
    return error_total;
}

amrex::Real normal_vector_neumann_test_impl(
    amr_wind::Field& vof, amr_wind::IntField& iblk_fld, const int& dir)
{
    const amrex::Real ref_m = 16.0_rt * (1.0_rt - 0.0_rt);
    amrex::Real error_total = 0.0_rt;

    for (int lev = 0; lev < vof.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            vof(lev), iblk_fld(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& vof_arr,
                amrex::Array4<int const> const& iblank) -> amrex::Real {
                amrex::Real error = 0.0_rt;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    int ibdy =
                        (iblank(i, j, k) != iblank(i - 1, j, k)) ? -1 : 0;
                    int jbdy =
                        (iblank(i, j, k) != iblank(i, j - 1, k)) ? -1 : 0;
                    int kbdy =
                        (iblank(i, j, k) != iblank(i, j, k - 1)) ? -1 : 0;
                    ibdy = (iblank(i, j, k) != iblank(i + 1, j, k)) ? +1 : ibdy;
                    jbdy = (iblank(i, j, k) != iblank(i, j + 1, k)) ? +1 : jbdy;
                    kbdy = (iblank(i, j, k) != iblank(i, j, k + 1)) ? +1 : kbdy;
                    amrex::Real mxn, myn, mzn;
                    amr_wind::multiphase::
                        youngs_finite_difference_normal_neumann(
                            i, j, k, ibdy, jbdy, kbdy, vof_arr, mxn, myn, mzn);

                    // Use L1 norm, check for 0
                    if (iblank(i, j, k) == 1) {
                        // Only interested in normals from field cells

                        // Slope in non-normal directions should be 0
                        if (dir != 0) {
                            error += std::abs(mxn - 0.0_rt);
                        }
                        if (dir != 1) {
                            error += std::abs(myn - 0.0_rt);
                        }
                        if (dir != 2) {
                            error += std::abs(mzn - 0.0_rt);
                        }
                        // Slope in normal direction, at center, should be same
                        if (dir == 0 && i == 2) {
                            error += std::abs(mxn - ref_m);
                        }
                        if (dir == 1 && j == 2) {
                            error += std::abs(myn - ref_m);
                        }
                        if (dir == 2 && k == 2) {
                            error += std::abs(mzn - ref_m);
                        }
                    }
                });

                return error;
            });
    }
    return error_total;
}

amrex::Real fit_plane_test_impl(amr_wind::Field& vof, const int dir)
{
    amrex::Real error_total = 0.0_rt;

    for (int lev = 0; lev < vof.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            vof(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& vof_arr)
                -> amrex::Real {
                amrex::Real error = 0.0_rt;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    int ii = (dir != 0 ? i : 0);
                    int jj = (dir != 1 ? j : 0);
                    int kk = (dir != 2 ? k : 0);
                    // Check multiphase cells
                    if (ii + jj + kk == 3) {
                        amrex::Real mx, my, mz, alpha;
                        amr_wind::multiphase::fit_plane(
                            i, j, k, vof_arr, mx, my, mz, alpha);

                        // Check slope
                        error += std::abs(mx - (dir != 0 ? 0.5_rt : 0.0_rt));
                        error += std::abs(my - (dir != 1 ? 0.5_rt : 0.0_rt));
                        error += std::abs(mz - (dir != 2 ? 0.5_rt : 0.0_rt));
                        // Check intercept
                        error += std::abs(alpha - 0.5_rt);
                    }
                });

                return error;
            });
    }
    return error_total;
}

amrex::Real fit_plane_test_impl_h(
    amr_wind::Field& vof, const amrex::Real vof_val, const int dir)
{
    amrex::Real error_total = 0.0_rt;
    const amrex::Real vv = vof_val;

    for (int lev = 0; lev < vof.repo().num_active_levels(); ++lev) {

        error_total += amrex::ReduceSum(
            vof(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& vof_arr)
                -> amrex::Real {
                amrex::Real error = 0.0_rt;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    int ii = (dir == 0 ? i : (dir == 1 ? j : k));
                    // Check multiphase cells
                    if (ii == 1) {
                        amrex::Real mx, my, mz, alpha;
                        amr_wind::multiphase::fit_plane(
                            i, j, k, vof_arr, mx, my, mz, alpha);

                        // Check slope
                        error += std::abs(mx - (dir == 0 ? 1.0_rt : 0.0_rt));
                        error += std::abs(my - (dir == 1 ? 1.0_rt : 0.0_rt));
                        error += std::abs(mz - (dir == 2 ? 1.0_rt : 0.0_rt));
                        // Check intercept
                        error += std::abs(alpha - vv);
                    }
                });

                return error;
            });
    }
    return error_total;
}

} // namespace

TEST_F(VOFOpTest, volume_intercept)
{
    // Initialize random number generator
    amrex::InitRandom(0);
    for (int n = 0; n < 100; ++n) {
        amrex::Real mx = amrex::Random();
        amrex::Real my = amrex::Random();
        amrex::Real mz = amrex::Random();
        amrex::Real vof = amrex::Random();
        // Scale slope values, like in fit_plane
        amrex::Real mm2 = mx + my + mz;
        mx = mx / mm2;
        my = my / mm2;
        mz = mz / mm2;
        // Limit vof values to multiphase range
        vof = amrex::max<amrex::Real>(
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt,
            amrex::min<amrex::Real>(
                1.0_rt -
                    (std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt),
                vof));
        // Get intercept value and check for nan
        amrex::Real alpha =
            amr_wind::multiphase::volume_intercept(mx, my, mz, vof);
        EXPECT_EQ(alpha, alpha);

        // Set one of the slope components to 0, then try again
        auto idx = (int)std::floor(amrex::Random() * 3.0_rt);
        switch (idx) {
        case 0:
            mx = 0.0_rt;
            break;
        case 1:
            my = 0.0_rt;
            break;
        default:
            mz = 0.0_rt;
            break;
        }
        // Scale slope values, like in fit_plane
        mm2 = mx + my + mz;
        mx = mx / mm2;
        my = my / mm2;
        mz = mz / mm2;
        // Get intercept value and check for nan
        alpha = amr_wind::multiphase::volume_intercept(mx, my, mz, vof);
        EXPECT_EQ(alpha, alpha);
    }

    // Check specific problem case
    amrex::Real alpha =
        amr_wind::multiphase::volume_intercept(0.5_rt, 0.5_rt, 0.0_rt, 0.5_rt);
    EXPECT_EQ(alpha, alpha);
}

TEST_F(VOFOpTest, interface_normal)
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
    const int ncomp = 1;
    const int nghost = 3;
    auto& vof = repo.declare_field("vof", ncomp, nghost);

    amrex::Real error_total = 0.0_rt;
    // constant in x
    init_vof(vof, 0);
    error_total = normal_vector_test_impl(vof, 0);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0_rt, tol);
    // constant in y
    init_vof(vof, 1);
    error_total = normal_vector_test_impl(vof, 1);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0_rt, tol);
    // constant in z
    init_vof(vof, 2);
    error_total = normal_vector_test_impl(vof, 2);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0_rt, tol);
}

TEST_F(VOFOpTest, interface_plane)
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
    const int ncomp = 1;
    const int nghost = 3;
    auto& vof = repo.declare_field("vof", ncomp, nghost);

    /* -- Diagonal plane, 2D orientation -- */
    amrex::Real error_total = 0.0_rt;
    // constant in x
    init_vof(vof, 0);
    error_total = fit_plane_test_impl(vof, 0);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0_rt, tol);
    // constant in y
    init_vof(vof, 1);
    error_total = fit_plane_test_impl(vof, 1);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0_rt, tol);
    // constant in z
    init_vof(vof, 2);
    error_total = fit_plane_test_impl(vof, 2);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0_rt, tol);

    /* -- "Horizontal" plane, 1D orientation -- */
    amrex::InitRandom(0);
    for (int n = 0; n < 20; ++n) {
        amrex::Real vof_val = amrex::Random();
        vof_val = amrex::max<amrex::Real>(
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt,
            amrex::min<amrex::Real>(
                1.0_rt -
                    (std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt),
                vof_val));
        // in x
        init_vof_h(vof, vof_val, 0);
        error_total = fit_plane_test_impl_h(vof, vof_val, 0);
        amrex::ParallelDescriptor::ReduceRealSum(error_total);
        EXPECT_NEAR(error_total, 0.0_rt, tol);
        // in y
        init_vof_h(vof, vof_val, 1);
        error_total = fit_plane_test_impl_h(vof, vof_val, 1);
        amrex::ParallelDescriptor::ReduceRealSum(error_total);
        EXPECT_NEAR(error_total, 0.0_rt, tol);
        // in z
        init_vof_h(vof, vof_val, 2);
        error_total = fit_plane_test_impl_h(vof, vof_val, 2);
        amrex::ParallelDescriptor::ReduceRealSum(error_total);
        EXPECT_NEAR(error_total, 0.0_rt, tol);
    }
}

TEST_F(VOFOpTest, interface_normal_neumann)
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
    const int ncomp = 1;
    const int nghost = 3;
    auto& vof = repo.declare_field("vof", ncomp, nghost);
    auto& iblank = repo.declare_int_field("iblank", ncomp, 1);

    // Check agreement / disagreement with ordinary calculation
    amrex::Real error_total = 0.0_rt;
    // iblank constant in x, vof varies in every direction
    init_vof(vof);
    init_iblank(iblank, 0);
    error_total = normal_vector_neumann_test_impl(vof, iblank);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0_rt, tol);
    // iblank constant in y, vof varies in every direction
    init_vof(vof);
    init_iblank(iblank, 1);
    error_total = normal_vector_neumann_test_impl(vof, iblank);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0_rt, tol);
    // iblank constant in z, vof varies in every direction
    init_vof(vof);
    init_iblank(iblank, 2);
    error_total = normal_vector_neumann_test_impl(vof, iblank);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0_rt, tol);

    // Confirm that neumann gets 0 when it should
    // iblank varies in y and z, vof varies in x
    init_vof_iblank(vof, iblank, 0);
    error_total = normal_vector_neumann_test_impl(vof, iblank, 0);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0_rt, tol);
    // iblank varies in x and z, vof varies in y
    init_vof_iblank(vof, iblank, 1);
    error_total = normal_vector_neumann_test_impl(vof, iblank, 1);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0_rt, tol);
    // iblank varies in x and y, vof varies in z
    init_vof_iblank(vof, iblank, 2);
    error_total = normal_vector_neumann_test_impl(vof, iblank, 2);
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    EXPECT_NEAR(error_total, 0.0_rt, tol);
}

} // namespace amr_wind_tests
