#include "abl_test_utils.H"
#include "src/wind_energy/MOData.H"
#include "src/wind_energy/ShearStress.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {
kynema_sgf::MOData make_mo(const amrex::Real alpha_h)
{
    kynema_sgf::MOData mo;
    mo.zref = 5.0_rt;
    mo.z0 = 0.1_rt;
    mo.z0t = 0.1_rt;
    mo.vmag_mean = 8.0_rt;
    mo.vel_mean[0] = 8.0_rt;
    mo.vel_mean[1] = 0.0_rt;
    mo.theta_mean = 300.0_rt;
    mo.alpha_h = alpha_h;
    return mo;
}

//! Temperature of the converged Monin-Obukhov profile at height z
amrex::Real mo_temperature(const kynema_sgf::MOData& mo, const amrex::Real z)
{
    return mo.surf_temp -
           (mo.surf_temp_flux * mo.phi_h(z) / (mo.utau * mo.kappa));
}

//! Run the heat-flux mode, then recover it from the near-surface temperature
void check_near_surface_round_trip(
    const amrex::Real heat_flux, const amrex::Real alpha_h)
{
    // update_fluxes stops once utau changes by less than 1e-5
    constexpr amrex::Real tol = 1.0e-3_rt;
    constexpr amrex::Real near_height = 2.0_rt;

    auto flux_mo = make_mo(alpha_h);
    flux_mo.alg_type = kynema_sgf::MOData::ThetaCalcType::HEAT_FLUX;
    flux_mo.surf_temp_flux = heat_flux;
    flux_mo.update_fluxes();

    auto near_mo = make_mo(alpha_h);
    near_mo.alg_type =
        kynema_sgf::MOData::ThetaCalcType::NEAR_SURFACE_TEMPERATURE;
    near_mo.near_surf_height = near_height;
    near_mo.near_surf_temp = mo_temperature(flux_mo, near_height);
    near_mo.update_fluxes();

    EXPECT_NEAR(near_mo.surf_temp_flux, heat_flux, tol * std::abs(heat_flux));
    EXPECT_NEAR(near_mo.surf_temp, flux_mo.surf_temp, tol);
    EXPECT_NEAR(near_mo.utau, flux_mo.utau, tol * flux_mo.utau);
    // The solution passes through both specified temperatures
    EXPECT_NEAR(
        mo_temperature(near_mo, near_height), near_mo.near_surf_temp, tol);
    EXPECT_NEAR(mo_temperature(near_mo, near_mo.zref), near_mo.theta_mean, tol);
}
} // namespace

TEST_F(ABLTest, mo_near_surface_temperature_unstable)
{
    check_near_surface_round_trip(0.05_rt, 1.0_rt);
}

TEST_F(ABLTest, mo_near_surface_temperature_stable)
{
    check_near_surface_round_trip(-0.01_rt, 1.0_rt);
}

TEST_F(ABLTest, mo_near_surface_temperature_alpha_h)
{
    check_near_surface_round_trip(0.05_rt, 0.74_rt);
}

TEST_F(ABLTest, mo_near_surface_temperature_neutral)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    auto mo = make_mo(1.0_rt);
    mo.alg_type = kynema_sgf::MOData::ThetaCalcType::NEAR_SURFACE_TEMPERATURE;
    mo.near_surf_height = 2.0_rt;
    mo.near_surf_temp = mo.theta_mean;
    mo.update_fluxes();
    EXPECT_NEAR(mo.surf_temp_flux, 0.0_rt, tol);
    EXPECT_NEAR(mo.surf_temp, mo.theta_mean, tol * mo.theta_mean);
}

TEST_F(ABLTest, mo_alpha_h_wall_model_flux)
{
    // The wall model must apply the heat flux the MO solve was given, for any
    // alpha_h
    constexpr amrex::Real tol = 1.0e-3_rt;
    constexpr amrex::Real heat_flux = 0.05_rt;
    for (const amrex::Real alpha_h : {1.0_rt, 0.74_rt, 1.35_rt}) {
        auto mo = make_mo(alpha_h);
        mo.alg_type = kynema_sgf::MOData::ThetaCalcType::HEAT_FLUX;
        mo.surf_temp_flux = heat_flux;
        mo.update_fluxes();

        const kynema_sgf::ShearStressMoeng moeng(mo);
        EXPECT_NEAR(
            -moeng.calc_theta(mo.vmag_mean, mo.theta_mean), heat_flux,
            tol * heat_flux);
    }
}

TEST_F(ABLTest, mo_alpha_h_surface_temperature)
{
    // alpha_h scales the temperature difference that carries a given flux
    constexpr amrex::Real tol = 1.0e-3_rt;
    auto mo1 = make_mo(1.0_rt);
    mo1.alg_type = kynema_sgf::MOData::ThetaCalcType::HEAT_FLUX;
    mo1.surf_temp_flux = 0.05_rt;
    mo1.update_fluxes();

    auto mo2 = make_mo(0.74_rt);
    mo2.alg_type = kynema_sgf::MOData::ThetaCalcType::SURFACE_TEMPERATURE;
    mo2.surf_temp =
        mo1.theta_mean + (0.74_rt * (mo1.surf_temp - mo1.theta_mean));
    mo2.update_fluxes();
    EXPECT_NEAR(mo2.surf_temp_flux, mo1.surf_temp_flux, tol * 0.05_rt);
}

} // namespace kynema_sgf_tests
