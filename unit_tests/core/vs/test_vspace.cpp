#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/core/vs/vector_space.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

namespace {
constexpr amrex::Real tol =
    std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

void test_vector_create_impl()
{
    namespace vs = amr_wind::vs;
    amrex::Gpu::DeviceScalar<amrex::Real> ds(0.0_rt);
    auto* ddata = ds.dataPtr();
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*unused*/) {
        auto gv1 = vs::Vector::ihat();
        auto gv2 = vs::Vector::jhat();
        auto gv3 = vs::Vector::khat();
        auto gv4 = (gv1 ^ gv2);

        ddata[0] = vs::mag((gv3 - gv4));
    });

    EXPECT_NEAR(ds.dataValue(), 0.0_rt, tol);
}

void test_tensor_create_impl()
{
    namespace vs = amr_wind::vs;
    amrex::Gpu::DeviceScalar<amrex::Real> ds(0.0_rt);
    auto* ddata = ds.dataPtr();
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*unused*/) {
        auto t1 = vs::yrot(90.0_rt);
        auto t2 = vs::zrot(90.0_rt);
        auto t3 = t2 & t1;
        auto qrot = vs::quaternion(vs::Vector::one(), 120.0_rt);

        ddata[0] = vs::mag((t3 - qrot));
    });

    EXPECT_NEAR(ds.dataValue(), 0.0_rt, tol);
}

void test_rotations_impl()
{
    namespace vs = amr_wind::vs;
    const amrex::Real angle = 45.0_rt;
    const auto ang = amr_wind::utils::radians(angle);
    const auto cval = std::cos(ang);
    const auto sval = std::sin(ang);

    auto xrot = vs::xrot(angle);
    auto yrot = vs::yrot(angle);
    auto zrot = vs::zrot(angle);

    auto ivec = vs::Vector::ihat();
    auto jvec = vs::Vector::jhat();
    auto kvec = vs::Vector::khat();

#define CHECK_ON_GPU(expr1, expr2)                                             \
    {                                                                          \
        amrex::Gpu::DeviceScalar<amrex::Real> ds(                              \
            std::numeric_limits<amrex::Real>::epsilon());                      \
        auto* dv = ds.dataPtr();                                               \
        amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int) {                      \
            auto v1 = expr1;                                                   \
            auto v2 = expr2;                                                   \
            dv[0] = vs::mag((v1 - v2));                                        \
        });                                                                    \
        EXPECT_NEAR(ds.dataValue(), 0.0_rt, tol)                               \
            << "LHS = " #expr1 "\nRHS = " #expr2 << '\n';                      \
    }

    CHECK_ON_GPU((xrot & ivec), ivec);
    CHECK_ON_GPU((xrot & jvec), (vs::Vector{0.0_rt, cval, -sval}));
    CHECK_ON_GPU((xrot & kvec), (vs::Vector{0.0_rt, sval, cval}));

    CHECK_ON_GPU((yrot & jvec), jvec);
    CHECK_ON_GPU((yrot & ivec), (vs::Vector{cval, 0.0_rt, sval}));
    CHECK_ON_GPU((yrot & kvec), (vs::Vector{-sval, 0.0_rt, cval}));

    CHECK_ON_GPU((zrot & kvec), kvec);
    CHECK_ON_GPU((zrot & ivec), (vs::Vector{cval, -sval, 0.0_rt}));
    CHECK_ON_GPU((zrot & jvec), (vs::Vector{sval, cval, 0.0_rt}));

    CHECK_ON_GPU((kvec & zrot), kvec);
    CHECK_ON_GPU((ivec & zrot), (vs::Vector{cval, sval, 0.0_rt}));
    CHECK_ON_GPU((jvec & zrot), (vs::Vector{-sval, cval, 0.0_rt}));

#undef CHECK_ON_GPU
}

void test_device_capture_impl()
{
    namespace vs = amr_wind::vs;
    auto v1 = vs::Vector::ihat();
    auto vexpected = vs::Vector::khat();
    amrex::Gpu::DeviceScalar<amrex::Real> ds(
        std::numeric_limits<amrex::Real>::epsilon());
    auto* dv = ds.dataPtr();

    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*unused*/) {
        auto v2 = vs::Vector::jhat();
        auto vout = v1 ^ v2;

        dv[0] = vs::mag((vout - vexpected));
    });
    EXPECT_NEAR(ds.dataValue(), 0.0_rt, tol);
}

void test_device_lists_impl()
{
    namespace vs = amr_wind::vs;
    amrex::Gpu::DeviceVector<vs::Vector> dvectors(3);
    auto* dv = dvectors.data();

    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int /*unused*/) {
        auto v1 = vs::Vector::ihat();
        auto v2 = vs::Vector::jhat();
        auto v3 = vs::Vector::khat();

        dv[0] = v2 ^ v3;
        dv[1] = v3 ^ v1;
        dv[2] = v1 ^ v2;
    });

    amrex::Vector<vs::Vector> hvectors(3);
    amrex::Vector<vs::Vector> htrue{
        vs::Vector::ihat(), vs::Vector::jhat(), vs::Vector::khat()};
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, dvectors.begin(), dvectors.end(),
        hvectors.begin());

    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(vs::mag(htrue[i] - hvectors[i]), 0.0_rt, tol);
    }
}

} // namespace

TEST(VectorSpace, vector_create)
{
    namespace vs = amr_wind::vs;

    auto v1 = vs::Vector::ihat();
    auto v2 = vs::Vector::jhat();
    auto v3 = vs::Vector::khat();
    auto v4 = v1 ^ v2;
    auto v5 = v1 ^ v3;
    auto v6 = v2 ^ v3;

    EXPECT_NEAR(v3.x(), v4.x(), tol);
    EXPECT_NEAR(v3.y(), v4.y(), tol);
    EXPECT_NEAR(v3.z(), v4.z(), tol);

    EXPECT_NEAR(-v2.x(), v5.x(), tol);
    EXPECT_NEAR(-v2.y(), v5.y(), tol);
    EXPECT_NEAR(-v2.z(), v5.z(), tol);

    EXPECT_NEAR(v1.x(), v6.x(), tol);
    EXPECT_NEAR(v1.y(), v6.y(), tol);
    EXPECT_NEAR(v1.z(), v6.z(), tol);

    test_vector_create_impl();
}

TEST(VectorSpace, tensor_create) { test_tensor_create_impl(); }

TEST(VectorSpace, vector_rotations) { test_rotations_impl(); }

TEST(VectorSpace, device_capture) { test_device_capture_impl(); }

TEST(VectorSpace, device_lists) { test_device_lists_impl(); }

TEST(VectorSpace, vector_ops)
{
    namespace vs = amr_wind::vs;
    const vs::Vector v1{10.0_rt, 20.0_rt, 0.0_rt};
    const vs::Vector v2{1.0_rt, 2.0_rt, 0.0_rt};
    const auto v21 = vs::mag_sqr(v2.unit());

    EXPECT_NEAR((v1 & v2), 50.0_rt * v21, tol);
    EXPECT_NEAR((v1 & vs::Vector::khat()), 0.0_rt, tol);
    EXPECT_NEAR(
        vs::mag_sqr((v1 + v2) - vs::Vector{11.0_rt, 22.0_rt, 0.0_rt}), 0.0_rt,
        tol);
}

} // namespace amr_wind_tests
