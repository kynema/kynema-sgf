#include "aw_test_utils/AmrexTest.H"
#include "amr-wind/incflo.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

namespace {

void init_density(amr_wind::Field& density, const int k_thresh = -3)
{
    const int nlevels = density.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& darrs = density(lev).arrays();

        amrex::ParallelFor(
            density(lev), density.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                darrs[nbx](i, j, k) = (k > k_thresh) ? k + 4 : 0.0_rt;
            });
    }
    amrex::Gpu::streamSynchronize();
}

amrex::Real get_fgz_sum(amr_wind::Field& src_term)
{
    amrex::Real Fgz_sum = 0.0_rt;

    for (int lev = 0; lev < src_term.repo().num_active_levels(); ++lev) {
        Fgz_sum += amrex::ReduceSum(
            src_term(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& Fg_arr) -> amrex::Real {
                amrex::Real Fgz_sum_fab = 0.0_rt;

                amrex::Loop(bx, [=, &Fgz_sum_fab](int i, int j, int k) {
                    Fgz_sum_fab += Fg_arr(i, j, k, 2);
                });

                return Fgz_sum_fab;
            });
    }
    amrex::ParallelDescriptor::ReduceRealSum(Fgz_sum);
    return Fgz_sum;
}

void fgtest_kernel(
    const amrex::Real Fz_ref,
    const int ncells,
    const int nz,
    const amr_wind::FieldState fstate,
    const bool make_ref_dens = false)
{
    incflo my_incflo;
    my_incflo.init_mesh();
    if (make_ref_dens) {
        // Create reference density field (only used if turned on via parser)
        auto& ref_dens =
            my_incflo.sim().repo().declare_field("reference_density", 1, 3, 1);
        init_density(ref_dens, nz / 2);
    }
    my_incflo.init_amr_wind_modules();
    auto& density = my_incflo.sim().repo().get_field("density").state(
        amr_wind::FieldState::NPH);
    auto& velocity = my_incflo.sim().repo().get_field("velocity");
    auto& grad_p = my_incflo.sim().repo().get_field("gp");
    auto& Fg_field = my_incflo.icns().fields().src_term;
    // Set density field
    init_density(density);
    // Set old density to unity, avoid NaN
    density.state(amr_wind::FieldState::Old).setVal(1.0_rt);
    // Set zero velocity
    velocity.setVal(0.0_rt);
    // Set zero pressure gradient
    grad_p.setVal(0.0_rt);

    // Calculate forcing terms
    my_incflo.icns().compute_source_term(fstate);

    // Check forcing average value
    const amrex::Real Fz_avg = get_fgz_sum(Fg_field) / ncells;
    EXPECT_NEAR(
        Fz_ref, Fz_avg, std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
}

} // namespace

class GravityForcingTest : public AmrexTest
{
protected:
    void populate_parameters()
    {
        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{m_nx, m_ny, m_nz}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", m_nx);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{{1.0_rt, 1.0_rt, 1.0_rt}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);

            amrex::Vector<int> periodic{{1, 1, 1}};
            pp.addarr("is_periodic", periodic);
        }
        {
            amrex::ParmParse pp("incflo");
            pp.add("use_godunov", 1);
        }
        {
            amrex::ParmParse pp("ICNS");
            amrex::Vector<std::string> srcstr{"GravityForcing"};
            pp.addarr("source_terms", srcstr);
        }
    }

    const amrex::Real m_rho_0 = 2.0_rt;
    const amrex::Real m_Fg = -9.81_rt;
    const int m_nx = 2;
    const int m_ny = 2;
    const int m_nz = 16;
};

TEST_F(GravityForcingTest, full_term_u)
{
    // High-level setup
    populate_parameters();
    {
        amrex::ParmParse pp("ICNS");
        pp.add("use_perturb_pressure", false);
    }
    // Modify gravity to make sure it works
    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<amrex::Real> grav{0.0_rt, 0.0_rt, -5.0_rt};
        pp.addarr("gravity", grav);
    }
    // Expected average gravity term
    amrex::Real Fg = -5.0_rt;
    // Test with ordinary gravity term (rho not included)
    fgtest_kernel(Fg, m_nx * m_ny * m_nz, m_nz, amr_wind::FieldState::Old);
}

TEST_F(GravityForcingTest, full_term_rhou)
{
    // High-level setup
    populate_parameters();
    {
        amrex::ParmParse pp("ICNS");
        pp.add("use_perturb_pressure", false);
    }
    // Modify gravity to make sure it works
    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<amrex::Real> grav{0.0_rt, 0.0_rt, -5.0_rt};
        pp.addarr("gravity", grav);
    }
    // Expected average gravity term
    amrex::Real Fg = -5.0_rt;
    amrex::Real fac = 0.0_rt;
    for (int k = 0; k < m_nz; ++k) {
        fac += k + 4;
    }
    Fg *= fac / m_nz;
    // Test with ordinary gravity term (rho included)
    fgtest_kernel(Fg, m_nx * m_ny * m_nz, m_nz, amr_wind::FieldState::New);
}

TEST_F(GravityForcingTest, perturb_const)
{
    const amrex::Real rho_ref = 0.5_rt;
    const amrex::Real gz = -9.81_rt;
    // High-level setup
    populate_parameters();
    {
        amrex::ParmParse pp("ICNS");
        pp.add("use_perturb_pressure", true);
    }
    // Modify gravity to make sure it works
    {
        amrex::ParmParse pp("incflo");
        pp.add("density", rho_ref);
    }
    // Expected average gravity term
    amrex::Real Fg = (1.0_rt - rho_ref) * gz / 1.0_rt;
    // Test with ordinary gravity term (rho not multiplied)
    fgtest_kernel(Fg, m_nx * m_ny * m_nz, m_nz, amr_wind::FieldState::Old);
}

TEST_F(GravityForcingTest, perturb_field)
{
    const amrex::Real gz = -9.81_rt;
    // High-level setup
    populate_parameters();
    {
        amrex::ParmParse pp("ICNS");
        pp.add("use_perturb_pressure", true);
    }
    // Expected average gravity term
    amrex::Real Fg = gz;
    amrex::Real fac = 0.0_rt;
    for (int k = 0; k < m_nz; ++k) {
        fac += (k > m_nz / 2.0_rt) ? 0.0_rt : k + 4.0_rt;
    }
    Fg *= fac / m_nz;
    // Test with ordinary gravity term (rho included)
    fgtest_kernel(
        Fg, m_nx * m_ny * m_nz, m_nz, amr_wind::FieldState::New, true);
}

} // namespace amr_wind_tests
