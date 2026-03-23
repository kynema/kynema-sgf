#include <string>
#include <numbers>
#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/wind_energy/actuator/aero/AirfoilTable.H"
#include "amr-wind/utilities/trig_ops.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {
namespace {

std::stringstream generate_txt_airfoil()
{
    using namespace ::amr_wind::utils;
    std::stringstream ss;
    ss << 6 << '\n';
    ss << -180.0_rt << " " << 0.0_rt << " " << 0.0_rt << " " << 0.0_rt << '\n';
    ss << -160.0_rt << " " << two_pi() * radians(20.0_rt) << " " << 0.0_rt
       << " " << 0.0_rt << '\n';
    ss << -20.0_rt << " " << -two_pi() * radians(20.0_rt) << " " << 0.0_rt
       << " " << 0.0_rt << '\n';
    ss << 20.0_rt << " " << two_pi() * radians(20.0_rt) << " " << 0.0_rt << " "
       << 0.0_rt << '\n';
    ss << 160.0_rt << " " << -two_pi() * radians(20.0_rt) << " " << 0.0_rt
       << " " << 0.0_rt << '\n';
    ss << 180.0_rt << " " << 0.0_rt << " " << 0.0_rt << " " << 0.0_rt << '\n';

    return ss;
}

std::stringstream generate_openfast_airfoil()
{
    std::stringstream ss;
    ss << "!........................................ " << '\n';
    ss << "! Table of aerodynamics coefficients " << '\n';
    ss << "        6   NumAlf            ! Number of data lines in the "
          "following table "
       << '\n';
    ss << "!    Alpha      Cl      Cd        Cm  " << '\n';
    ss << "!    (deg)      (-)     (-)       (-) " << '\n';
    ss << "   -180.00    0.000   0.0407   0.0000 " << '\n';
    ss << "   -175.00    0.223   0.0507   0.0937 " << '\n';
    ss << "   -170.00    0.405   0.1055   0.1702 " << '\n';
    ss << "   -160.00    0.658   0.2982   0.2819 " << '\n';
    ss << "   -155.00    0.733   0.4121   0.3213 " << '\n';
    ss << "   -150.00    0.778   0.5308   0.3520 " << '\n';
    ss << "   -145.00    0.795   0.6503   0.3754 " << '\n';
    ss << "   -140.00    0.787   0.7672   0.3926 " << '\n';
    ss << "   -135.00    0.757   0.8785   0.4046 " << '\n';

    return ss;
}

} // namespace

TEST(Airfoil, read_txt_file)
{
    using AirfoilLoader = ::amr_wind::actuator::AirfoilLoader;
    auto ss = generate_txt_airfoil();

    auto af = AirfoilLoader::load_text_file(ss);
    EXPECT_EQ(af->num_entries(), 6);
    EXPECT_NEAR(
        af->aoa().front(), -1.0_rt * std::numbers::pi_v<amrex::Real>,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
    EXPECT_NEAR(
        af->aoa().back(), 1.0_rt * std::numbers::pi_v<amrex::Real>,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
}

TEST(Airfoil, read_openfast_file)
{
    using AirfoilLoader = ::amr_wind::actuator::AirfoilLoader;
    auto ss = generate_openfast_airfoil();

    auto af = AirfoilLoader::load_openfast_airfoil(ss);
    EXPECT_EQ(af->num_entries(), 6);
    EXPECT_NEAR(
        af->aoa().front(), -1.0_rt * std::numbers::pi_v<amrex::Real>,
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
    EXPECT_NEAR(
        af->aoa().back(), ::amr_wind::utils::radians(-150.0_rt),
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
}

TEST(Airfoil, airfoil_lookup)
{
    using AirfoilLoader = ::amr_wind::actuator::AirfoilLoader;
    auto ss = generate_txt_airfoil();

    auto af = AirfoilLoader::load_text_file(ss);
    EXPECT_EQ(af->num_entries(), 6);

    amrex::Real cl, cd;
    amrex::Vector<amrex::Real> aoa_test{-15.0_rt, -5.0_rt, 0.0_rt,
                                        5.0_rt,   10.0_rt, 15.0_rt};
    for (const auto& aoa_deg : aoa_test) {
        const amrex::Real aoa_rad = ::amr_wind::utils::radians(aoa_deg);
        (*af)(aoa_rad, cl, cd);
        EXPECT_NEAR(cl, ::amr_wind::utils::two_pi() * aoa_rad, 1.0e-3_rt);
        EXPECT_NEAR(cd, 0.0_rt, 1.0e-3_rt);
    }
}

} // namespace amr_wind_tests
