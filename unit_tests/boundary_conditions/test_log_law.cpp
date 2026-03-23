#include "gtest/gtest.h"
#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/boundary_conditions/wall_models/LogLaw.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

class LogLawTest : public MeshTest
{
protected:
    void populate_parameters() override { MeshTest::populate_parameters(); }
    [[nodiscard]] amrex::Real log_law_actual(const amrex::Real utau) const
    {
        return utau * (std::log(m_zref * utau / m_nu) / 0.384_rt + 4.27_rt);
    }

    const amrex::Real m_zref = 1.0_rt / 32.0_rt;
    // Molecular viscosity
    const amrex::Real m_nu = 1.0e-4_rt;
    // Set test tolerance to convergence test in LogLaw
    const amrex::Real m_tol = 1.0e-5_rt;
};

TEST_F(LogLawTest, test_log_law)
{
    populate_parameters();
    amr_wind::LogLaw ll;
    ll.zref = m_zref;
    ll.nu = m_nu;

    amrex::Real utau_expected = 0.1_rt;
    amrex::Real wspd = log_law_actual(utau_expected);
    auto utau = ll.get_utau(wspd);
    EXPECT_NEAR(utau, 0.1_rt, m_tol);

    utau_expected = 0.05_rt;
    EXPECT_NEAR(
        ll.get_utau(log_law_actual(utau_expected)), utau_expected, m_tol);

    utau_expected = 2.0_rt;
    EXPECT_NEAR(
        ll.get_utau(log_law_actual(utau_expected)), utau_expected, m_tol);

    utau_expected = 0.5_rt;
    ll.wspd_mean = log_law_actual(utau_expected);
    ll.update_utau_mean();
    EXPECT_NEAR(ll.utau_mean, utau_expected, m_tol);
}
} // namespace amr_wind_tests
