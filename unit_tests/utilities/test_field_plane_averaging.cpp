#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "AMReX_Box.H"
#include "AMReX_BoxArray.H"
#include "AMReX_BoxList.H"
#include "AMReX_Geometry.H"
#include "AMReX_RealBox.H"
#include "AMReX_Vector.H"
#include "amr-wind/utilities/FieldPlaneAveraging.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

class FieldPlaneAveragingTest : public MeshTest
{
public:
    void test_dir(int /*dir*/);
};

TEST_F(FieldPlaneAveragingTest, test_constant)
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

    // test the average of a constant is the same constant
    for (int dir = 0; dir < 3; ++dir) {

        amr_wind::FieldPlaneAveraging pa(velocityf, sim().time(), dir);
        pa();

        amrex::Real z = 0.5_rt * (problo[dir] + probhi[dir]);

        amrex::Real u = pa.line_average_interpolated(z, 0);
        amrex::Real v = pa.line_average_interpolated(z, 1);
        amrex::Real w = pa.line_average_interpolated(z, 2);

        EXPECT_NEAR(u0, u, tol);
        EXPECT_NEAR(v0, v, tol);
        EXPECT_NEAR(w0, w, tol);

        for (int i = 0; i < 8; ++i) {
            EXPECT_NEAR(0.0_rt, pa.line_derivative_of_average_cell(i, 0), tol);
            EXPECT_NEAR(0.0_rt, pa.line_derivative_of_average_cell(i, 1), tol);
            EXPECT_NEAR(0.0_rt, pa.line_derivative_of_average_cell(i, 2), tol);
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
    // const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> xlo =
    // geom.ProbLoArray(); const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx
    // = geom.CellSizeArray();
    const auto xlo = geom.ProbLoArray();
    const auto dx = geom.CellSizeArray();

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

TEST_F(FieldPlaneAveragingTest, test_linear)
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

    amr_wind::FieldPlaneAveraging pa(velocityf, sim().time(), dir);
    pa();

    constexpr int n = 20;
    const amrex::Real L = probhi[dir] - problo[dir];
    const amrex::Real dx = L / ((amrex::Real)n);
    const amrex::Real hchLo =
        problo[dir] + (0.5_rt * mesh().Geom(0).CellSizeArray()[dir]);
    const amrex::Real hchHi =
        probhi[dir] - (0.5_rt * mesh().Geom(0).CellSizeArray()[dir]);

    // test along a line at n equidistant points
    for (int i = 0; i < n; ++i) {
        const amrex::Real x = problo[dir] + (i * dx);

        const amrex::Array<amrex::Real, 3> u = {
            pa.line_average_interpolated(x, 0),
            pa.line_average_interpolated(x, 1),
            pa.line_average_interpolated(x, 2)};

        amrex::Real xtest;

        if (x < hchLo) {
            // test near bottom boundary where solution is not linearly
            // interpolated but instead copied from first cell
            xtest = hchLo;
        } else if (x > hchHi) {
            // test near top boundary where solution is not linearly
            // interpolated but instead copied from last cell
            xtest = hchHi;
        } else {
            // interior is linearly interpolated
            xtest = x;
        }

        // test each velocity field u = u0 + u0*x
        for (int j = 0; j < 3; ++j) {
            EXPECT_NEAR(u0[j] * (xtest + 1.0_rt), u[j], tol);
        }
    }

    for (int i = 0; i < 8; ++i) {
        EXPECT_NEAR(u0[0], pa.line_derivative_of_average_cell(i, 0), tol);
        EXPECT_NEAR(u0[1], pa.line_derivative_of_average_cell(i, 1), tol);
        EXPECT_NEAR(u0[2], pa.line_derivative_of_average_cell(i, 2), tol);
    }
}

namespace {

void add_periodic(
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

        for (int d = 0; d < 3; ++d) {
            if (d != dir) {
                velocity(i, j, k, 0) += std::cos(a[d] * x[d]);
                velocity(i, j, k, 1) += std::sin(a[d] * x[d]);
                velocity(i, j, k, 2) +=
                    std::sin(a[d] * x[d]) * std::cos(a[d] * x[d]);
            }
        }
    });
}

} // namespace

void FieldPlaneAveragingTest::test_dir(int dir)
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
    a[0] = periods * amr_wind::utils::two_pi() / (probhi[0] - problo[0]);
    a[1] = periods * amr_wind::utils::two_pi() / (probhi[1] - problo[1]);
    a[2] = periods * amr_wind::utils::two_pi() / (probhi[2] - problo[2]);

    run_algorithm(
        mesh().num_levels(), velocity,
        [&](const int lev, const amrex::MFIter& mfi) {
            auto vel = velocity[lev]->array(mfi);
            const auto& bx = mfi.validbox();

            add_periodic(dir, a, mesh().Geom(lev), bx, vel);
        });

    amr_wind::FieldPlaneAveraging pa(velocityf, sim().time(), dir);
    pa();

    amrex::Real x = 0.5_rt * (problo[dir] + probhi[dir]);
    amrex::Real u = pa.line_average_interpolated(x, 0);
    amrex::Real v = pa.line_average_interpolated(x, 1);
    amrex::Real w = pa.line_average_interpolated(x, 2);

    // test that a periodic function orthogonal to dir
    // averages out to a constant
    EXPECT_NEAR(u0, u, tol);
    EXPECT_NEAR(v0, v, tol);
    EXPECT_NEAR(w0, w, tol);
}

TEST_F(FieldPlaneAveragingTest, test_xdir) { test_dir(0); }
TEST_F(FieldPlaneAveragingTest, test_ydir) { test_dir(1); }
TEST_F(FieldPlaneAveragingTest, test_zdir) { test_dir(2); }

} // namespace amr_wind_tests
