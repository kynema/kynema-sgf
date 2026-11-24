#include <array>
#include <numbers>
#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/iter_tools.H"
#include "ks_test_utils/test_utils.H"
#include "src/ocean_waves/relaxation_zones/stokes_waves_K.H"
#include "src/utilities/math_ops.H"
#include "src/utilities/constants.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

class WaveTheoriesTest : public MeshTest
{};

TEST_F(WaveTheoriesTest, StokesWaves)
{
    constexpr amrex::Real coeff_tol = 1.0e-4_rt;

    constexpr amrex::Real wavenumber = 2.0_rt;
    constexpr amrex::Real water_depth = 0.376991_rt;
    constexpr int stokes_order = 5;
    amrex::Real c0, a11, a22, b22, c2, a31, a33, b31, a42;
    amrex::Real a44, b42, b44, c4, a51, a53, a55, b53, b55;

    kynema_sgf::ocean_waves::relaxation_zones::stokes_coefficients(
        stokes_order, wavenumber, water_depth, c0, a11, a22, b22, c2, a31, a33,
        b31, a42, a44, b42, b44, c4, a51, a53, a55, b53, b55);

    // Gold coefficient values taken from table 2 of
    // Fenton, J. Fifth Order Stokes Theory for Steady Waves
    // Journal of Waterway, Port, Coastal and Ocean Engineering, 1985, 111,
    // 216-234

    const amrex::Real gold_A11 = 1.208490_rt;
    const amrex::Real gold_A22 = 0.799840_rt;
    const amrex::Real gold_A31 = -9.105340_rt;
    const amrex::Real gold_A33 = 0.368275_rt;
    const amrex::Real gold_A42 = -12.196150_rt;
    const amrex::Real gold_A44 = 0.058723_rt;
    const amrex::Real gold_A51 = 108.46831725_rt;
    const amrex::Real gold_A53 = -6.941756_rt;
    const amrex::Real gold_A55 = -0.074979_rt;
    const amrex::Real gold_B22 = 2.502414_rt;
    const amrex::Real gold_B31 = -5.731666_rt;
    const amrex::Real gold_B42 = -32.407508_rt;
    const amrex::Real gold_B44 = 14.033758_rt;
    const amrex::Real gold_B53 = -103.44536875_rt;
    const amrex::Real gold_B55 = 37.200027_rt;
    const amrex::Real gold_C0 = 0.798448_rt;
    const amrex::Real gold_C2 = 1.940215_rt;
    const amrex::Real gold_C4 = -12.970403_rt;

    EXPECT_NEAR(gold_A11, a11, coeff_tol);
    EXPECT_NEAR(gold_A22, a22, coeff_tol);
    EXPECT_NEAR(gold_A31, a31, coeff_tol);
    EXPECT_NEAR(gold_A33, a33, coeff_tol);
    EXPECT_NEAR(gold_A42, a42, coeff_tol);
    EXPECT_NEAR(gold_A44, a44, coeff_tol);
    EXPECT_NEAR(gold_A51, a51, coeff_tol);
    EXPECT_NEAR(gold_A53, a53, coeff_tol);
    EXPECT_NEAR(gold_A55, a55, coeff_tol);
    EXPECT_NEAR(gold_B22, b22, coeff_tol);
    EXPECT_NEAR(gold_B31, b31, coeff_tol);
    EXPECT_NEAR(gold_B42, b42, coeff_tol);
    EXPECT_NEAR(gold_B44, b44, coeff_tol);
    EXPECT_NEAR(gold_B53, b53, coeff_tol);
    EXPECT_NEAR(gold_B55, b55, coeff_tol);
    EXPECT_NEAR(gold_C0, c0, coeff_tol);
    EXPECT_NEAR(gold_C2, c2, coeff_tol);
    EXPECT_NEAR(gold_C4, c4, coeff_tol);
}

TEST_F(WaveTheoriesTest, StokesWavesFreeSurfaceProfile)
{
    constexpr amrex::Real tol = 1.0e-4_rt;

    constexpr amrex::Real g = 9.81_rt;
    constexpr int stokes_order = 5;
    // wavenumber k and water_depth d chosen so that kd = 0.758932_rt
    // to match value of column 3 from table 2 in
    // Fenton, J. Fifth Order Stokes Theory for Steady Waves
    // Journal of Waterway, Port, Coastal and Ocean Engineering, 1985, 111,
    // 216-234
    amrex::Real wavenumber = 2.0_rt;
    amrex::Real water_depth = 0.376991_rt;
    amrex::Real wavelength =
        2.0_rt * std::numbers::pi_v<amrex::Real> / wavenumber;
    amrex::Real wave_height = 0.1_rt;
    amrex::Real zsl = 0.0_rt;
    amrex::Real x = 0.0_rt;
    amrex::Real z = -0.25_rt;
    amrex::Real phase_offset = 0.0_rt;
    amrex::Real time = 0.0_rt;
    amrex::Real eta = 0.0_rt;
    amrex::Real u_w = 0.0_rt;
    amrex::Real v_w = 0.0_rt;
    amrex::Real w_w = 0.0_rt;

    kynema_sgf::ocean_waves::relaxation_zones::stokes_waves(
        stokes_order, wavelength, water_depth, wave_height, zsl, g, x, z, time,
        phase_offset, kynema_sgf::constants::LOW_NUM,
        kynema_sgf::constants::LARGE_NUM, eta, u_w, v_w, w_w);

    // Coefficients values taken from column 3 of table 2 of
    // Fenton, J. Fifth Order Stokes Theory for Steady Waves
    // Journal of Waterway, Port, Coastal and Ocean Engineering, 1985, 111,
    // 216-234
    amrex::Real B22 = 2.502414_rt;
    amrex::Real B31 = -5.731666_rt;
    amrex::Real B42 = -32.407508_rt;
    amrex::Real B44 = 14.033758_rt;
    amrex::Real B53 = -103.44536875_rt;
    amrex::Real B55 = 37.200027_rt;

    amrex::Real eps = wavenumber * wave_height / 2.0_rt;
    amrex::Real S = 2.0_rt * std::exp(2.0_rt * wavenumber * water_depth) /
                    (std::exp(4.0_rt * wavenumber * water_depth) + 1.0_rt);
    amrex::Real C = 1.0_rt - S;
    amrex::Real C0 = std::sqrt(std::tanh(wavenumber * water_depth));
    amrex::Real C2 = C0 * (2.0_rt + (7.0_rt * kynema_sgf::utils::powi(S, 2))) /
                     (4.0_rt * kynema_sgf::utils::powi(C, 2));
    amrex::Real C4 =
        C0 *
        (4.0_rt + (32.0_rt * S) - (116.0_rt * kynema_sgf::utils::powi(S, 2)) -
         (400.0_rt * kynema_sgf::utils::powi(S, 3)) -
         (71.0_rt * kynema_sgf::utils::powi(S, 4)) +
         (146.0_rt * kynema_sgf::utils::powi(S, 5))) /
        (32.0_rt * kynema_sgf::utils::powi(C, 5));
    amrex::Real wave_speed = (C0 + (kynema_sgf::utils::powi(eps, 2) * C2) +
                              (kynema_sgf::utils::powi(eps, 4) * C4)) *
                             std::sqrt(g / wavenumber);

    amrex::Real omega = wave_speed * wavenumber;
    amrex::Real phase = (wavenumber * x) - (omega * time) - phase_offset;

    // Check against Eq. (14) from Fenton 1985
    amrex::Real eta_theory =
        (((eps * std::cos(phase)) +
          (kynema_sgf::utils::powi(eps, 2) * B22 * std::cos(2.0_rt * phase)) +
          (kynema_sgf::utils::powi(eps, 3) * B31 *
           (std::cos(phase) - std::cos(3.0_rt * phase))) +
          (kynema_sgf::utils::powi(eps, 4) *
           ((B42 * std::cos(2.0_rt * phase)) +
            (B44 * std::cos(4.0_rt * phase)))) +
          (kynema_sgf::utils::powi(eps, 5) *
           ((-(B53 + B55) * std::cos(phase)) +
            (B53 * std::cos(3.0_rt * phase)) +
            (B55 * std::cos(5.0_rt * phase))))) /
         wavenumber) +
        zsl;

    EXPECT_NEAR(eta, eta_theory, tol);

    // Re-evaluate with new set of coefficients
    // Deep-water limit (k*d->\infty)
    water_depth = 100.0_rt;
    wave_height = 0.16_rt;
    wavenumber = 0.156_rt;
    wavelength = 2.0_rt * std::numbers::pi_v<amrex::Real> / wavenumber;
    zsl = 0.0_rt;
    x = 4.0_rt;
    z = 0.0_rt;
    phase_offset = std::numbers::pi_v<amrex::Real>;
    time = 2.7_rt;
    eta = 0.0_rt;
    u_w = 0.0_rt;
    v_w = 0.0_rt;
    w_w = 0.0_rt;

    kynema_sgf::ocean_waves::relaxation_zones::stokes_waves(
        stokes_order, wavelength, water_depth, wave_height, zsl, g, x, z, time,
        phase_offset, kynema_sgf::constants::LOW_NUM,
        kynema_sgf::constants::LARGE_NUM, eta, u_w, v_w, w_w);

    // Coefficients values taken from column 1 of table 2 of
    // Fenton, J. Fifth Order Stokes Theory for Steady Waves
    // Journal of Waterway, Port, Coastal and Ocean Engineering, 1985, 111,
    // 216-234
    B22 = 0.5_rt;
    B31 = -0.375_rt;
    B42 = 0.3333333_rt;
    B44 = 0.3333333_rt;
    B53 = 0.7734375_rt;
    B55 = 0.3255208_rt;

    eps = wavenumber * wave_height / 2.0_rt;
    // Coefficients computed analytically by taking limit kd -> \infty and
    // simplified accordingly
    // Note that in this limit S = 0 and C = 1 and thus they are omitted here
    C0 = 1.0_rt;
    C2 = 0.5_rt;
    C4 = 0.125_rt;
    wave_speed = (C0 + (kynema_sgf::utils::powi(eps, 2) * C2) +
                  (kynema_sgf::utils::powi(eps, 4) * C4)) *
                 std::sqrt(g / wavenumber);

    omega = wave_speed * wavenumber;
    phase = (wavenumber * x) - (omega * time) - phase_offset;

    // Matches Eq. (18) from Fenton 1985
    eta_theory =
        (((eps * std::cos(phase)) +
          (kynema_sgf::utils::powi(eps, 2) * B22 * std::cos(2.0_rt * phase)) +
          (kynema_sgf::utils::powi(eps, 3) * B31 *
           (std::cos(phase) - std::cos(3.0_rt * phase))) +
          (kynema_sgf::utils::powi(eps, 4) *
           ((B42 * std::cos(2.0_rt * phase)) +
            (B44 * std::cos(4.0_rt * phase)))) +
          (kynema_sgf::utils::powi(eps, 5) *
           ((-(B53 + B55) * std::cos(phase)) +
            (B53 * std::cos(3.0_rt * phase)) +
            (B55 * std::cos(5.0_rt * phase))))) /
         wavenumber) +
        zsl;

    EXPECT_NEAR(eta, eta_theory, tol);
}

TEST_F(WaveTheoriesTest, StokesWavesVelocityComponents)
{
    constexpr amrex::Real tol = 1.0e-4_rt;

    constexpr amrex::Real g = 9.81_rt;
    constexpr int stokes_order = 5;
    // wavenumber k and water_depth d chosen so that kd = 0.758932_rt
    // to match value of column 3 from table 2 in
    // Fenton, J. Fifth Order Stokes Theory for Steady Waves
    // Journal of Waterway, Port, Coastal and Ocean Engineering, 1985, 111,
    // 216-234
    constexpr amrex::Real wavenumber = 2.0_rt;
    constexpr amrex::Real water_depth = 0.376991_rt;
    constexpr amrex::Real wavelength =
        2.0_rt * std::numbers::pi_v<amrex::Real> / wavenumber;
    constexpr amrex::Real wave_height = 0.1_rt;
    constexpr amrex::Real zsl = 0.0_rt;
    constexpr amrex::Real x = 0.0_rt;
    constexpr amrex::Real z = -0.25_rt;
    constexpr amrex::Real phase_offset = 0.0_rt;
    constexpr amrex::Real time = 0.0_rt;
    amrex::Real eta = 0.0_rt;
    amrex::Real u_w = 0.0_rt;
    amrex::Real v_w = 0.0_rt;
    amrex::Real w_w = 0.0_rt;

    kynema_sgf::ocean_waves::relaxation_zones::stokes_waves(
        stokes_order, wavelength, water_depth, wave_height, zsl, g, x, z, time,
        phase_offset, kynema_sgf::constants::LOW_NUM,
        kynema_sgf::constants::LARGE_NUM, eta, u_w, v_w, w_w);

    // Coefficients values taken from column 3 of table 2 of
    // Fenton, J. Fifth Order Stokes Theory for Steady Waves
    // Journal of Waterway, Port, Coastal and Ocean Engineering, 1985, 111,
    // 216-234
    const amrex::Real A11 = 1.208490_rt;
    const amrex::Real A22 = 0.799840_rt;
    const amrex::Real A31 = -9.105340_rt;
    const amrex::Real A33 = 0.368275_rt;
    const amrex::Real A42 = -12.196150_rt;
    const amrex::Real A44 = 0.058723_rt;
    const amrex::Real A51 = 108.46831725_rt;
    const amrex::Real A53 = -6.941756_rt;
    const amrex::Real A55 = -0.074979_rt;

    const amrex::Real eps = wavenumber * wave_height / 2.0_rt;
    const amrex::Real S =
        2.0_rt * std::exp(2.0_rt * wavenumber * water_depth) /
        (std::exp(4.0_rt * wavenumber * water_depth) + 1.0_rt);
    const amrex::Real C = 1.0_rt - S;
    const amrex::Real C0 = std::sqrt(std::tanh(wavenumber * water_depth));
    const amrex::Real C2 = C0 *
                           (2.0_rt + (7.0_rt * kynema_sgf::utils::powi(S, 2))) /
                           (4.0_rt * kynema_sgf::utils::powi(C, 2));
    const amrex::Real C4 =
        C0 *
        (4.0_rt + (32.0_rt * S) - (116.0_rt * kynema_sgf::utils::powi(S, 2)) -
         (400.0_rt * kynema_sgf::utils::powi(S, 3)) -
         (71.0_rt * kynema_sgf::utils::powi(S, 4)) +
         (146.0_rt * kynema_sgf::utils::powi(S, 5))) /
        (32.0_rt * kynema_sgf::utils::powi(C, 5));
    const amrex::Real wave_speed =
        (C0 + (kynema_sgf::utils::powi(eps, 2) * C2) +
         (kynema_sgf::utils::powi(eps, 4) * C4)) *
        std::sqrt(g / wavenumber);

    const amrex::Real omega = wave_speed * wavenumber;
    const amrex::Real phase = (wavenumber * x) - (omega * time) - phase_offset;

    // Compare with theoretical Results from Kinnas
    // https://www.sciencedirect.com/science/article/pii/S0029801817306066
    // Define coefficients using Eq.(19)
    amrex::Vector<amrex::Real> a(stokes_order);
    a[0] = A11 + ((eps * eps) * A31) + (kynema_sgf::utils::powi(eps, 4) * A51);
    a[1] = A22 + ((eps * eps) * A42);
    a[2] = A33 + ((eps * eps) * A53);
    a[3] = A44;
    a[4] = A55;

    // Horizontal velocity from Eq.(21) and vertical velocity from Eq.(23) in
    // Kinnas
    // Initialize to first order terms
    amrex::Real horizontal_velocity =
        eps * a[0] * std::cosh(wavenumber * (water_depth + (z - zsl))) *
        std::cos(phase);
    amrex::Real vertical_velocity =
        eps * a[0] * std::sinh(wavenumber * (water_depth + (z - zsl))) *
        std::sin(phase);

    for (int n = 1; n < stokes_order; ++n) {
        horizontal_velocity +=
            kynema_sgf::utils::powi(eps, (n + 1)) * (n + 1) * a[n] *
            std::cosh((n + 1) * wavenumber * (water_depth + (z - zsl))) *
            std::cos((n + 1) * phase);
        vertical_velocity +=
            kynema_sgf::utils::powi(eps, (n + 1)) * (n + 1) * a[n] *
            std::sinh((n + 1) * wavenumber * (water_depth + (z - zsl))) *
            std::sin((n + 1) * phase);
    }
    horizontal_velocity *= (C0 * std::sqrt(g / wavenumber));
    vertical_velocity *= (C0 * std::sqrt(g / wavenumber));

    EXPECT_NEAR(u_w, horizontal_velocity, tol);
    EXPECT_NEAR(w_w, vertical_velocity, tol);
}

TEST_F(WaveTheoriesTest, StokesWaveLength)
{
    // Values of wave_height, wave_period and water_depth taken from
    // https://www.sciencedirect.com/science/article/pii/S0029801817306066
    amrex::Real wave_height = 0.16_rt;
    amrex::Real wave_period = 1.6_rt;
    amrex::Real water_depth = 18.0_rt;
    int wave_order = 2;
    constexpr amrex::Real g = 9.81_rt;
    constexpr amrex::Real tol_lambda =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt;
    // To return computed wavelength from first guess of wavenumber k
    int iter_max = -1;

    // Check initial guess of wavenumber k for Newton iterations
    amrex::Real lambda =
        kynema_sgf::ocean_waves::relaxation_zones::stokes_wave_length(
            wave_period, water_depth, wave_height, wave_order, g, tol_lambda,
            iter_max);

    const amrex::Real k_newton =
        2.0_rt * std::numbers::pi_v<amrex::Real> / lambda;

    // Compare with expected wavenumber from theory k = omega^2/g, where omega =
    // 2Pi/wave_period
    const amrex::Real k_theory =
        (2.0_rt * std::numbers::pi_v<amrex::Real> / wave_period) *
        (2.0_rt * std::numbers::pi_v<amrex::Real> / wave_period) / g;

    EXPECT_NEAR(
        k_newton, k_theory, std::numeric_limits<float>::epsilon() * 1.0e1_rt);

    // Check wave theory
    wave_height = 0.2_rt;
    wave_period = 2.2_rt;
    water_depth = 5.0_rt;
    wave_order = 5;
    iter_max = 20;

    lambda = kynema_sgf::ocean_waves::relaxation_zones::stokes_wave_length(
        wave_period, water_depth, wave_height, wave_order, g, tol_lambda,
        iter_max);

    // Relation to check is Eq.(24) from course notes:
    // https://www.caee.utexas.edu/prof/kinnas/ce358/oenotes/kinnas_stokes11.pdf
    amrex::Real k = 2.0_rt * std::numbers::pi_v<amrex::Real> / lambda;
    const amrex::Real RHS1 = 2.0_rt * std::numbers::pi_v<amrex::Real> /
                             (wave_period * std::sqrt(g * k));

    amrex::Real S = 1.0_rt / std::cosh(2.0_rt * k * water_depth);
    amrex::Real C = 1.0_rt - S;
    amrex::Real eps = k * wave_height / 2.0_rt;

    amrex::Real C0 = std::sqrt(std::tanh(k * water_depth));
    amrex::Real C2 = C0 * (2.0_rt + (7.0_rt * S * S)) / (4.0_rt * C * C);
    const amrex::Real C4 =
        C0 *
        (4.0_rt + (32.0_rt * S) - (116.0_rt * kynema_sgf::utils::powi(S, 2)) -
         (400.0_rt * kynema_sgf::utils::powi(S, 3)) -
         (71.0_rt * kynema_sgf::utils::powi(S, 4)) +
         (146.0_rt * kynema_sgf::utils::powi(S, 5))) /
        (32.0_rt * kynema_sgf::utils::powi(C, 5));
    const amrex::Real LHS1 = C0 + (kynema_sgf::utils::powi(eps, 2) * C2) +
                             (kynema_sgf::utils::powi(eps, 4) * C4);
    EXPECT_NEAR(
        LHS1, RHS1, std::numeric_limits<amrex::Real>::epsilon() * 1.0e8_rt);

    // Reevaluate with a new set of conditions
    wave_height = 0.05_rt;
    wave_period = 1.5_rt;
    water_depth = 0.9_rt;
    wave_order = 3;
    iter_max = 31;
    lambda = kynema_sgf::ocean_waves::relaxation_zones::stokes_wave_length(
        wave_period, water_depth, wave_height, wave_order, g, tol_lambda,
        iter_max);

    k = 2.0_rt * std::numbers::pi_v<amrex::Real> / lambda;
    const amrex::Real RHS2 = 2.0_rt * std::numbers::pi_v<amrex::Real> /
                             (wave_period * std::sqrt(g * k));

    S = 1.0_rt / std::cosh(2.0_rt * k * water_depth);
    C = 1.0_rt - S;
    eps = k * wave_height / 2.0_rt;

    C0 = std::sqrt(std::tanh(k * water_depth));
    C2 = C0 * (2.0_rt + (7.0_rt * S * S)) / (4.0_rt * C * C);
    const amrex::Real LHS2 = C0 + (kynema_sgf::utils::powi(eps, 2) * C2);
    EXPECT_NEAR(
        LHS2, RHS2, std::numeric_limits<amrex::Real>::epsilon() * 1.0e8_rt);
}

namespace {
//! Period-averaged volume flux per unit span at x = 0, integrated from the bed
//! to the free surface, for Stokes waves with the given Eulerian mean current
amrex::Real stokes_mean_volume_flux(
    const int order,
    const amrex::Real wavelength,
    const amrex::Real depth,
    const amrex::Real height,
    const amrex::Real g,
    const amrex::Real mean_current)
{
    namespace rz = kynema_sgf::ocean_waves::relaxation_zones;
    const amrex::Real c = rz::stokes_wave_speed(
        order, wavelength, depth, height, g, mean_current);
    const amrex::Real period = wavelength / c;
    constexpr int nt = 400;
    constexpr int nz = 4000;
    amrex::Real flux = 0.0_rt;
    for (int it = 0; it < nt; ++it) {
        const amrex::Real t = period * it / nt;
        amrex::Real eta{0.0_rt}, u{0.0_rt}, v{0.0_rt}, w{0.0_rt};
        rz::stokes_waves(
            order, wavelength, depth, height, 0.0_rt, g, 0.0_rt, 0.0_rt, t,
            0.0_rt, kynema_sgf::constants::LOW_NUM,
            kynema_sgf::constants::LARGE_NUM, eta, u, v, w, mean_current);
        // Midpoint rule from the bed (z = -depth) to the surface eta
        const amrex::Real dz = (eta + depth) / nz;
        amrex::Real col = 0.0_rt;
        for (int iz = 0; iz < nz; ++iz) {
            const amrex::Real z = -depth + ((iz + 0.5_rt) * dz);
            amrex::Real eta_z{0.0_rt};
            rz::stokes_waves(
                order, wavelength, depth, height, 0.0_rt, g, 0.0_rt, z, t,
                0.0_rt, kynema_sgf::constants::LOW_NUM,
                kynema_sgf::constants::LARGE_NUM, eta_z, u, v, w, mean_current);
            col += u * dz;
        }
        flux += col / nt;
    }
    return flux;
}
} // namespace

TEST_F(WaveTheoriesTest, StokesZeroMassTransportLinearLimit)
{
    // At second order the return current is the linear mass transport divided
    // by the depth: -(eps^2 / 2) sqrt(coth(kd)) sqrt(g / k^3) / d
    namespace rz = kynema_sgf::ocean_waves::relaxation_zones;
    constexpr amrex::Real g = 9.81_rt;
    for (const amrex::Real depth : {0.5_rt, 2.0_rt, 4.0_rt}) {
        const amrex::Real wavelength = 4.0_rt;
        const amrex::Real height = 0.1_rt;
        const amrex::Real k =
            2.0_rt * std::numbers::pi_v<amrex::Real> / wavelength;
        const amrex::Real eps = k * height / 2.0_rt;
        const amrex::Real expected = -0.5_rt * eps * eps *
                                     std::sqrt(1.0_rt / std::tanh(k * depth)) *
                                     std::sqrt(g / (k * k * k)) / depth;
        EXPECT_NEAR(
            rz::stokes_zero_mass_transport_current(
                2, wavelength, depth, height, g),
            expected, 1.0e-12_rt);
        EXPECT_LT(expected, 0.0_rt);
    }
}

TEST_F(WaveTheoriesTest, StokesZeroMassTransportFlux)
{
    // Without the return current the waves carry a net volume flux equal to
    // the Stokes transport -d u1; with it, the net flux vanishes up to the
    // truncation of the fifth-order theory
    namespace rz = kynema_sgf::ocean_waves::relaxation_zones;
    constexpr amrex::Real g = 9.81_rt;
    constexpr int order = 5;
    // (wavelength, depth, height): steep deep-water wave as in the regular wave
    // cases, and a moderate intermediate-depth wave
    const std::array<std::array<amrex::Real, 3>, 2> cases{
        {{4.0_rt, 4.0_rt, 0.3585_rt}, {4.0_rt, 1.0_rt, 0.15_rt}}};
    for (const auto& cs : cases) {
        const amrex::Real u1 = rz::stokes_zero_mass_transport_current(
            order, cs[0], cs[1], cs[2], g);
        const amrex::Real q = -u1 * cs[1];
        const amrex::Real flux0 =
            stokes_mean_volume_flux(order, cs[0], cs[1], cs[2], g, 0.0_rt);
        const amrex::Real flux1 =
            stokes_mean_volume_flux(order, cs[0], cs[1], cs[2], g, u1);
        EXPECT_GT(q, 0.0_rt);
        EXPECT_NEAR(flux0, q, 0.03_rt * q);
        EXPECT_NEAR(flux1, 0.0_rt, 0.03_rt * q);
    }
}

TEST_F(WaveTheoriesTest, StokesWaveLengthZeroMassTransport)
{
    // Wavelength from the period without mass transport satisfies the
    // dispersion relation with the return current, and is shorter than the
    // wavelength for zero mean current (the waves travel against the current)
    namespace rz = kynema_sgf::ocean_waves::relaxation_zones;
    constexpr amrex::Real g = 9.81_rt;
    constexpr amrex::Real tol = 1.0e-10_rt;
    constexpr int order = 5;
    const amrex::Real period = 1.6_rt;
    const amrex::Real depth = 4.0_rt;
    const amrex::Real height = 0.3585_rt;
    const amrex::Real L0 =
        rz::stokes_wave_length(period, depth, height, order, g, tol, 40);
    const amrex::Real L1 = rz::stokes_wave_length_zero_mass_transport(
        period, depth, height, order, g, tol, 40);
    ASSERT_GT(L0, 0.0_rt);
    ASSERT_GT(L1, 0.0_rt);
    EXPECT_LT(L1, L0);
    const amrex::Real u1 =
        rz::stokes_zero_mass_transport_current(order, L1, depth, height, g);
    const amrex::Real c1 =
        rz::stokes_wave_speed(order, L1, depth, height, g, u1);
    EXPECT_NEAR(L1 / c1, period, 1.0e-8_rt);
    // Zero-current wave speed from the same function matches the existing
    // dispersion solution
    const amrex::Real c0 =
        rz::stokes_wave_speed(order, L0, depth, height, g, 0.0_rt);
    EXPECT_NEAR(L0 / c0, period, 1.0e-8_rt);
}

} // namespace kynema_sgf_tests
