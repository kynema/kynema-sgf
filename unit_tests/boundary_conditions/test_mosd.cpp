#include "gtest/gtest.h"
#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/boundary_conditions/wall_models/MOSD.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

class MOSDTest : public MeshTest
{
protected:
    void populate_parameters() override { MeshTest::populate_parameters(); }

    [[nodiscard]] amrex::Real get_u_loglaw(const amrex::Real utau) const
    {
        return utau *
               (1.0_rt / 0.4_rt * std::log(m_zref * utau / m_nu) + 5.0_rt);
    }

    // vertical location (z = dx)
    const amrex::Real m_zref = 6.2831853_rt / 128.0_rt;
    // Molecular viscosity
    const amrex::Real m_nu = 1.0e-5_rt;
    // Set test tolerance to convergence test in MOSD
    const amrex::Real m_tol = 1.0e-5_rt;
};

TEST_F(MOSDTest, test_mosd)
{
    populate_parameters();
    amr_wind::MOSD md;

    md.amplitude = 0.05_rt;
    md.wavenumber = 4;
    md.omega = 0.8_rt;
    md.time = 0;

    const amrex::Real utau = 0.1_rt;
    const amrex::Real u_dx = get_u_loglaw(utau);
    const amrex::Real v_dx = 0.0_rt;
    const amrex::Real unit_nor = 0.0_rt;
    {
        const amrex::Real tau_wave_expected = 0.0_rt;
        const amrex::Real x_c = std::numbers::pi_v<amrex::Real>;
        EXPECT_NEAR(
            md.get_dyn_tau(u_dx, v_dx, x_c, unit_nor), tau_wave_expected,
            m_tol);
    }

    {
        const amrex::Real tau_wave_expected = 0.0_rt;
        const amrex::Real x_c = 0.392699_rt; // This give sin(pi/two) = one//
        // However, because this location is at the negative slope of the wave,
        // the heaviside function will make the tau=zero
        EXPECT_NEAR(
            md.get_dyn_tau(u_dx, v_dx, x_c, unit_nor), tau_wave_expected,
            m_tol);
    }

    {
        const amrex::Real tau_wave_expected = 0.02493_rt;
        const amrex::Real x_c = 1.0_rt; // This give sin(four)= negative value
        // At this location we are in the positive slope, therefore heaviside
        // function will allow for stress
        EXPECT_NEAR(
            md.get_dyn_tau(u_dx, v_dx, x_c, unit_nor), tau_wave_expected,
            m_tol);
    }
}
} // namespace amr_wind_tests
