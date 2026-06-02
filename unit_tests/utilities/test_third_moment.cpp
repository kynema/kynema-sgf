#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/iter_tools.H"

#include "AMReX_Box.H"
#include "AMReX_BoxArray.H"
#include "AMReX_BoxList.H"
#include "AMReX_Geometry.H"
#include "AMReX_RealBox.H"
#include "AMReX_Vector.H"

#include "src/utilities/FieldPlaneAveraging.H"
#include "src/utilities/ThirdMomentAveraging.H"
#include "src/utilities/trig_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

class ThirdMomentAveragingTest : public MeshTest
{
public:
    void test_dir(int /*dir*/);
};

TEST_F(ThirdMomentAveragingTest, test_constant)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    constexpr amrex::Real u0 = 2.3_rt, v0 = 3.5_rt, w0 = 5.6_rt;

    populate_parameters();
    initialize_mesh();

    auto& frepo = mesh().field_repo();
    auto& velocityf = frepo.declare_field("velocity", 3);

    auto velocity = velocityf.vec_ptrs();

    // initialize level 0 to a constant
    velocity[0]->setVal(u0, 0, 1);
    velocity[0]->setVal(v0, 1, 1);
    velocity[0]->setVal(w0, 2, 1);

    const auto& problo = mesh().Geom(0).ProbLoArray();
    const auto& probhi = mesh().Geom(0).ProbHiArray();

    // test the fluctuation third-moments of a constant are zero
    for (int dir = 0; dir < 3; ++dir) {

        kynema_sgf::FieldPlaneAveraging pa(velocityf, sim().time(), dir, 0);
        pa();

        kynema_sgf::ThirdMomentAveraging uuu(pa, pa, pa);
        uuu();

        amrex::Real x = 0.5_rt * (problo[dir] + probhi[dir]);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    EXPECT_NEAR(
                        uuu.line_average_interpolated(x, i, j, k), 0.0_rt, tol);
                }
            }
        }
    }
}

namespace {

void add_linear(
    int dir,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> a,
    const amrex::Geometry& geom,
    const amrex::Box& bx,
    const amrex::Array4<amrex::Real>& velocity)
{
    auto xlo = geom.ProbLoArray();
    auto dx = geom.CellSizeArray();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::GpuArray<amrex::Real, 3> x = {
            xlo[0] + ((i + 0.5_rt) * dx[0]), xlo[1] + ((j + 0.5_rt) * dx[1]),
            xlo[2] + ((k + 0.5_rt) * dx[2])};
        velocity(i, j, k, 0) += a[0] * x[dir];
        velocity(i, j, k, 1) += a[1] * x[dir];
        velocity(i, j, k, 2) += a[2] * x[dir];
    });
}

} // namespace

TEST_F(ThirdMomentAveragingTest, test_linear)
{

    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> u0 = {
        {1.0_rt, 3.5_rt, 5.6_rt}};

    populate_parameters();
    initialize_mesh();

    auto& frepo = mesh().field_repo();
    auto& velocityf = frepo.declare_field("velocity", 3);
    auto velocity = velocityf.vec_ptrs();

    velocity[0]->setVal(u0[0], 0, 1);
    velocity[0]->setVal(u0[1], 1, 1);
    velocity[0]->setVal(u0[2], 2, 1);

    constexpr int dir = 2;

    const auto& problo = mesh().Geom(0).ProbLoArray();
    const auto& probhi = mesh().Geom(0).ProbHiArray();

    run_algorithm(
        mesh().num_levels(), velocity,
        [&](const int lev, const amrex::MFIter& mfi) {
            auto vel = velocity[lev]->array(mfi);
            const auto& bx = mfi.validbox();
            add_linear(dir, u0, mesh().Geom(0), bx, vel);
        });

    kynema_sgf::FieldPlaneAveraging pa(velocityf, sim().time(), dir, 0);
    pa();
    kynema_sgf::ThirdMomentAveraging uuu(pa, pa, pa);
    uuu();

    constexpr int n = 20;
    const amrex::Real L = probhi[dir] - problo[dir];
    const amrex::Real dx = L / ((amrex::Real)n);

    // test along a line at n equidistant points
    for (int i = 0; i < n; ++i) {

        const amrex::Real x = problo[dir] + (i * dx);

        for (int m = 0; m < 3; ++m) {
            for (int ncomp = 0; ncomp < 3; ++ncomp) {
                for (int p = 0; p < 3; ++p) {
                    EXPECT_NEAR(
                        0.0_rt, uuu.line_average_interpolated(x, m, ncomp, p),
                        tol);
                }
            }
        }
    }
}

namespace {

void add_periodic(
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> a,
    const amrex::Geometry& geom,
    const amrex::Box& bx,
    const amrex::Array4<amrex::Real>& velocity)
{
    auto xlo = geom.ProbLoArray();
    auto dx = geom.CellSizeArray();

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::GpuArray<amrex::Real, 3> x = {
            xlo[0] + ((i + 0.5_rt) * dx[0]), xlo[1] + ((j + 0.5_rt) * dx[1]),
            xlo[2] + ((k + 0.5_rt) * dx[2])};
        for (int d = 0; d < 3; ++d) {
            velocity(i, j, k, 0) += std::cos(a[d] * x[d]);
            velocity(i, j, k, 1) += std::sin(a[d] * x[d]);
            velocity(i, j, k, 2) +=
                std::sin(a[d] * x[d]) * std::cos(a[d] * x[d]);
        }
    });
}

} // namespace

void ThirdMomentAveragingTest::test_dir(int dir)
{

    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    constexpr amrex::Real u0 = 2.3_rt, v0 = 3.5_rt, w0 = 5.6_rt;

    populate_parameters();
    initialize_mesh();

    auto& frepo = mesh().field_repo();
    auto& velocityf = frepo.declare_field("velocity", 3);
    auto velocity = velocityf.vec_ptrs();

    velocity[0]->setVal(u0, 0, 1);
    velocity[0]->setVal(v0, 1, 1);
    velocity[0]->setVal(w0, 2, 1);

    constexpr int periods = 3;

    const auto& problo = mesh().Geom(0).ProbLoArray();
    const auto& probhi = mesh().Geom(0).ProbHiArray();

    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> a;
    a[0] = periods * kynema_sgf::utils::two_pi() / (probhi[0] - problo[0]);
    a[1] = periods * kynema_sgf::utils::two_pi() / (probhi[1] - problo[1]);
    a[2] = periods * kynema_sgf::utils::two_pi() / (probhi[2] - problo[2]);

    run_algorithm(
        mesh().num_levels(), velocity,
        [&](const int lev, const amrex::MFIter& mfi) {
            auto vel = velocity[lev]->array(mfi);
            const auto& bx = mfi.validbox();

            add_periodic(a, mesh().Geom(lev), bx, vel);
        });

    kynema_sgf::FieldPlaneAveraging pa(velocityf, sim().time(), dir, 0);
    pa();

    kynema_sgf::ThirdMomentAveraging uuu(pa, pa, pa);
    uuu();

    amrex::Real x =
        (0.5_rt + 0.01_rt * amrex::Random()) * (problo[dir] + probhi[dir]);

    // Non-zero third moments are permutations of <u' v' w'> = 1/4.
    EXPECT_NEAR(0.25_rt, uuu.line_average_interpolated(x, 0, 1, 2), tol);
    EXPECT_NEAR(0.25_rt, uuu.line_average_interpolated(x, 0, 2, 1), tol);
    EXPECT_NEAR(0.25_rt, uuu.line_average_interpolated(x, 1, 0, 2), tol);
    EXPECT_NEAR(0.25_rt, uuu.line_average_interpolated(x, 1, 2, 0), tol);
    EXPECT_NEAR(0.25_rt, uuu.line_average_interpolated(x, 2, 0, 1), tol);
    EXPECT_NEAR(0.25_rt, uuu.line_average_interpolated(x, 2, 1, 0), tol);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                const bool is_012_perm =
                    (i != j) && (j != k) && (i != k) && (i + j + k == 3);
                if (!is_012_perm) {
                    EXPECT_NEAR(
                        0.0_rt, uuu.line_average_interpolated(x, i, j, k), tol);
                }
            }
        }
    }
}

TEST_F(ThirdMomentAveragingTest, test_xdir) { test_dir(0); }
TEST_F(ThirdMomentAveragingTest, test_ydir) { test_dir(1); }
TEST_F(ThirdMomentAveragingTest, test_zdir) { test_dir(2); }

} // namespace kynema_sgf_tests
