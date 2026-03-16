#include "gtest/gtest.h"
#include "aw_test_utils/MeshTest.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include "aw_test_utils/test_utils.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

namespace {

amrex::Real get_val_at_kindex(
    amr_wind::Field& field,
    amr_wind::Field& divisor,
    const int comp,
    const int kref)
{
    const int lev = 0;
    amrex::Real error_total = 0;

    error_total += amrex::ReduceSum(
        field(lev), divisor(lev), 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx, amrex::Array4<amrex::Real const> const& f_arr,
            amrex::Array4<amrex::Real const> const& div_arr) -> amrex::Real {
            amrex::Real error = 0;

            amrex::Loop(bx, [=, &error](int i, int j, int k) {
                // Check if current cell is just above lower wall
                if (k == kref) {
                    // Add field value to output
                    error += std::sqrt(f_arr(i, j, k, comp) / div_arr(i, j, k));
                }
            });

            return error;
        });
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    return error_total;
}

void init_field3(amr_wind::Field& fld, amrex::Real srate)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();

    amrex::Real offset = 0.0_rt;
    if (fld.field_location() == amr_wind::FieldLoc::CELL) {
        offset = 0.5_rt;
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = fld(lev).arrays();

        amrex::ParallelFor(
            fld(lev), fld.num_grow(), fld.num_comp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                const amrex::Real z = problo[2] + ((k + offset) * dx[2]);
                farrs[nbx](i, j, k, n) = (z / 2.0_rt * srate) + 2.0_rt;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void init_field1(amr_wind::Field& fld, amrex::Real tgrad)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();

    amrex::Real offset = 0.0_rt;
    if (fld.field_location() == amr_wind::FieldLoc::CELL) {
        offset = 0.5_rt;
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = fld(lev).arrays();

        amrex::ParallelFor(
            fld(lev), fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real z = problo[2] + ((k + offset) * dx[2]);

                farrs[nbx](i, j, k, 0) = z * tgrad;
            });
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace

class TurbLESTestBC : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{10, 20, 30}};
            pp.addarr("n_cell", ncell);
            pp.add("blocking_factor", 2);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{{10.0_rt, 10.0_rt, 10.0_rt}};
            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
            amrex::Vector<int> periodic{{1, 1, 0}};
            pp.addarr("is_periodic", periodic);
            // Boundary conditions
            amrex::ParmParse ppzhi("zhi");
            ppzhi.add("type", (std::string) "slip_wall");
            // zlo is defined in each case
        }
    }
    void one_eq_ksgs_setup_params() const
    {
        {
            amrex::ParmParse pp("turbulence");
            pp.add("model", (std::string) "OneEqKsgsM84");
        }
        {
            amrex::ParmParse pp("OneEqKsgsM84_coeffs");
            pp.add("Ceps", m_Ceps);
            pp.add("Ce", m_Ce);
        }
        {
            amrex::ParmParse pp("incflo");
            amrex::Vector<std::string> physics{"ABL"};
            pp.addarr("physics", physics);
            pp.add("density", m_rho0);
            amrex::Vector<amrex::Real> vvec{8.0_rt, 0.0_rt, 0.0_rt};
            pp.addarr("velocity", vvec);
            amrex::Vector<amrex::Real> gvec{0.0_rt, 0.0_rt, -m_gravz};
            pp.addarr("gravity", gvec);
        }
        {
            amrex::ParmParse pp("ABL");
            amrex::Vector<amrex::Real> t_hts{0.0_rt, 100.0_rt, 400.0_rt};
            pp.addarr("temperature_heights", t_hts);
            amrex::Vector<amrex::Real> t_vals{265.0_rt, 265.0_rt, 268.0_rt};
            pp.addarr("temperature_values", t_vals);
        }
        {
            amrex::ParmParse pp("transport");
            pp.add("viscosity", m_mu);
            pp.add("reference_temperature", m_Tref);
        }
    }
    void test_calls_body(bool do_postsolve = false)
    {
        // Initialize necessary parts of solver
        initialize_mesh();
        auto& pde_mgr = sim().pde_manager();
        auto& icns_eq = pde_mgr.register_icns();
        sim().init_physics();
        sim().create_turbulence_model();
        icns_eq.initialize();

        // Get turbulence model
        auto& tmodel = sim().turbulence_model();
        // Set up velocity field with constant strainrate
        auto& vel = sim().repo().get_field("velocity");
        init_field3(vel, m_srate);
        // Set up uniform unity density field
        auto& dens = sim().repo().get_field("density");
        dens.setVal(m_rho0);
        // Set up temperature field with constant gradient in z
        auto& temp = sim().repo().get_field("temperature");
        init_field1(temp, m_Tgz);
        // Give values to tlscale and tke arrays
        auto& tlscale = sim().repo().get_field("turb_lscale");
        tlscale.setVal(m_tlscale_val);
        auto& tke = sim().repo().get_field("tke");
        tke.setVal(m_tke_val);

        // Perform post init for physics
        // (turns on wall model if specified)
        for (auto& pp : sim().physics()) {
            pp->post_init_actions();
        }
        // Perform fillpatch to allow BCs to operate
        vel.fillpatch(0.0_rt);
        temp.fillpatch(0.0_rt);
        // Ensure nonzero viscosity (only molecular at this point)
        icns_eq.compute_mueff(amr_wind::FieldState::New);

        // Advance states (important for wall model)
        pde_mgr.advance_states();
        // Compute diffusion term to use BCs (wall model)
        icns_eq.compute_diffusion_term(amr_wind::FieldState::New);

        // Post_solve afterward
        if (do_postsolve) {
            icns_eq.post_solve_actions();
        }

        // Update turbulent viscosity directly
        tmodel.update_turbulent_viscosity(
            amr_wind::FieldState::New, DiffusionType::Crank_Nicolson);
    }

    const amrex::Real m_dx = 10.0_rt / 10.0_rt;
    const amrex::Real m_dy = 10.0_rt / 20.0_rt;
    const amrex::Real m_dz = 10.0_rt / 30.0_rt;
    const amrex::Real m_tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    // Parser inputs for turbulence model
    const amrex::Real m_Ceps = 0.11_rt;
    const amrex::Real m_Ce = 0.99_rt;
    const amrex::Real m_Tref = 263.5_rt;
    const amrex::Real m_gravz = 10.0_rt;
    const amrex::Real m_rho0 = 1.2_rt;
    // Constants for fields
    const amrex::Real m_srate = 20.0_rt;
    const amrex::Real m_Tgz = 2.0_rt;
    const amrex::Real m_tlscale_val = 1.1_rt;
    const amrex::Real m_tke_val = 0.1_rt;
    // Molecular viscosity
    const amrex::Real m_mu = 1.0e-2_rt;
};

TEST_F(TurbLESTestBC, test_1eqKsgs_noslip)
{
    populate_parameters();
    {
        amrex::ParmParse pp("zlo");
        pp.add("type", (std::string) "no_slip_wall");
    }
    one_eq_ksgs_setup_params();
    test_calls_body();
    auto& muturb = sim().repo().get_field("mu_turb");

    // Get shear production field
    auto& shear_prod = sim().repo().get_field("shear_prod");

    // Check for constant value in bulk
    auto shear_bulk =
        get_val_at_kindex(shear_prod, muturb, 0, 1) / 10.0_rt / 20.0_rt;
    EXPECT_NEAR(shear_bulk, m_srate, m_tol);

    // Velocity gradients assumed in setup (init_field3)
    const amrex::Real uz_bulk = m_srate / 2.0_rt;

    // Naive cell-centered answer, with no_slip_wall (Dirichlet)
    const amrex::Real uz_wallcell =
        (uz_bulk * 1.5_rt * m_dz + 2.0_rt - 0.0_rt) / (2.0_rt * m_dz);
    // Wall-normal direction
    const amrex::Real wz_wallcell = uz_wallcell;

    const amrex::Real s_naive = std::sqrt(
        (2.0_rt * wz_wallcell * wz_wallcell) +
        (2.0_rt * uz_wallcell * uz_wallcell));
    // Check for different value just above wall due to BC
    auto shear_wall =
        get_val_at_kindex(shear_prod, muturb, 0, 0) / 10.0_rt / 20.0_rt;
    // Check that the result is not equal to the naive value
    EXPECT_GT(std::abs(shear_wall - s_naive), m_tol);
    // Answer that accounts for location of wall at cell face
    const amrex::Real uz_wallface =
        ((uz_bulk * 1.5_rt * m_dz + 2.0_rt) / 3.0_rt + uz_bulk * 0.5_rt * m_dz +
         2.0_rt - 0.0_rt) /
        m_dz;
    const amrex::Real wz_wallface = uz_wallface;
    const amrex::Real s_true = std::sqrt(
        (2.0_rt * wz_wallface * wz_wallface) +
        (2.0_rt * uz_wallface * uz_wallface));
    EXPECT_NEAR(shear_wall, s_true, m_tol);
}

TEST_F(TurbLESTestBC, test_1eqKsgs_slip)
{
    populate_parameters();
    {
        amrex::ParmParse pp("zlo");
        pp.add("type", (std::string) "slip_wall");
    }
    one_eq_ksgs_setup_params();
    test_calls_body();
    auto& muturb = sim().repo().get_field("mu_turb");

    // Get shear production field
    auto& shear_prod = sim().repo().get_field("shear_prod");

    // Check for constant value in bulk
    auto shear_bulk =
        get_val_at_kindex(shear_prod, muturb, 0, 1) / 10.0_rt / 20.0_rt;
    EXPECT_NEAR(shear_bulk, m_srate, m_tol);

    // Check for different value just above wall due to BC
    auto shear_wall =
        get_val_at_kindex(shear_prod, muturb, 0, 0) / 10.0_rt / 20.0_rt;
    // Answer that accounts for slip_wall BC
    // (tangential extrapolation, dirichlet normal)
    const amrex::Real uz_bulk = m_srate / 2.0_rt;
    const amrex::Real wz_wall = ((uz_bulk * 1.5_rt * m_dz + 2.0_rt) / 3.0_rt +
                                 uz_bulk * 0.5_rt * m_dz + 2.0_rt - 0.0_rt) /
                                m_dz;
    const amrex::Real s_true =
        std::sqrt((2.0_rt * wz_wall * wz_wall) + (2.0_rt * uz_bulk * uz_bulk));
    EXPECT_NEAR(shear_wall, s_true, m_tol);
}

TEST_F(TurbLESTestBC, test_1eqKsgs_wallmodel)
{
    constexpr amrex::Real kappa = 0.4_rt;
    constexpr amrex::Real z0 = 0.11_rt;
    populate_parameters();
    {
        amrex::ParmParse pp("zlo");
        pp.add("type", (std::string) "wall_model");
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("wall_shear_stress_type", (std::string) "local");
        pp.add("kappa", kappa);
        pp.add("surface_roughness_z0", z0);
    }
    {
        // Explicit diffusion populates wall cells, doesn't update velocity
        amrex::ParmParse pp("incflo");
        pp.add("diffusion_type", 0);
    }
    one_eq_ksgs_setup_params();
    const bool do_postsolve = true;
    test_calls_body(do_postsolve);
    auto& muturb = sim().repo().get_field("mu_turb");

    // Get shear production field
    auto& shear_prod = sim().repo().get_field("shear_prod");

    // Check for constant value in bulk
    auto shear_bulk =
        get_val_at_kindex(shear_prod, muturb, 0, 1) / 10.0_rt / 20.0_rt;
    EXPECT_NEAR(shear_bulk, m_srate, m_tol);

    // Velocity gradients assumed in setup (init_field3)
    const amrex::Real uz_bulk = m_srate / 2.0_rt;
    // Get value just above wall due to BC
    auto shear_wall =
        get_val_at_kindex(shear_prod, muturb, 0, 0) / 10.0_rt / 20.0_rt;
    // (tangential extrapolation, dirichlet normal)
    const amrex::Real wz_wall = ((uz_bulk * 1.5_rt * m_dz + 2.0_rt) / 3.0_rt +
                                 uz_bulk * 0.5_rt * m_dz + 2.0_rt - 0.0_rt) /
                                m_dz;
    const amrex::Real s_true =
        std::sqrt((2.0_rt * wz_wall * wz_wall) + (2.0_rt * uz_bulk * uz_bulk));
    EXPECT_NEAR(shear_wall, s_true, m_tol);
}

TEST_F(TurbLESTestBC, test_1eqKsgs_wallmodel_failnofillpatch)
{
    constexpr amrex::Real kappa = 0.4_rt;
    constexpr amrex::Real z0 = 0.11_rt;
    populate_parameters();
    {
        amrex::ParmParse pp("zlo");
        pp.add("type", (std::string) "wall_model");
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("wall_shear_stress_type", (std::string) "local");
        pp.add("kappa", kappa);
        pp.add("aerodynamic_roughness_length", z0);
        pp.add("thermal_roughness_length", z0 + 0.02_rt);
    }
    {
        // Explicit diffusion populates wall cells, doesn't update velocity
        amrex::ParmParse pp("incflo");
        pp.add("diffusion_type", 0);
    }
    one_eq_ksgs_setup_params();
    const bool do_postsolve = false;
    test_calls_body(do_postsolve);
    auto& muturb = sim().repo().get_field("mu_turb");

    // Get shear production field
    auto& shear_prod = sim().repo().get_field("shear_prod");

    // Check for constant value in bulk
    auto shear_bulk =
        get_val_at_kindex(shear_prod, muturb, 0, 1) / 10.0_rt / 20.0_rt;
    EXPECT_NEAR(shear_bulk, m_srate, m_tol);

    // Velocity gradients assumed in setup (init_field3)
    const amrex::Real uz_bulk = m_srate / 2.0_rt;

    // Get value just above wall due to BC
    auto shear_wall =
        get_val_at_kindex(shear_prod, muturb, 0, 0) / 10.0_rt / 20.0_rt;
    // (tangential extrapolation, dirichlet normal)
    const amrex::Real wz_wall = ((uz_bulk * 1.5_rt * m_dz + 2.0_rt) / 3.0_rt +
                                 uz_bulk * 0.5_rt * m_dz + 2.0_rt - 0.0_rt) /
                                m_dz;
    const amrex::Real s_true =
        std::sqrt((2.0_rt * wz_wall * wz_wall) + (2.0_rt * uz_bulk * uz_bulk));
    // This is checking the correct value -- without the fillpatch, it is wrong
    // (passing test indicates fillpatch is needed after diffusion performed)
    EXPECT_GT(std::abs(shear_wall - s_true), m_tol);

    // Shear value at the wall used by wall model
    const amrex::Real zref = 0.5_rt * m_dz;
    const amrex::Real uref = (zref * uz_bulk) + 2.0_rt;
    const amrex::Real vmag_ref = std::sqrt(2.0_rt * uref * uref);
    const amrex::Real utau = kappa * vmag_ref / (std::log(zref / z0));
    const amrex::Real uz_wm =
        uref / vmag_ref * std::pow(utau, 2.0_rt) * m_rho0 / m_mu;

    // Velocity gradient with wallmodel value included as dirichlet
    const amrex::Real uz_wmdirichlet =
        (1.0_rt / 3.0_rt * (uz_bulk * 1.5_rt * m_dz + 2.0_rt) +
         1.0_rt * (uz_bulk * 0.5_rt * m_dz + 2.0_rt) -
         4.0_rt / 3.0_rt * uz_wm) /
        m_dz;

    // Naive answer, assumes wall Dirichlet
    const amrex::Real s_naive = std::sqrt(
        (2.0_rt * wz_wall * wz_wall) +
        (2.0_rt * uz_wmdirichlet * uz_wmdirichlet));
    // This is checking the wrong value
    // (passing test indicates it is wrong in expected way w/out fillpatch)
    EXPECT_NEAR(shear_wall, s_naive, std::abs(shear_wall) * m_tol);
}

TEST_F(TurbLESTestBC, test_1eqKsgs_zerogradient)
{
    populate_parameters();
    {
        amrex::ParmParse pp("zlo");
        pp.add("type", (std::string) "zero_gradient");
    }
    one_eq_ksgs_setup_params();
    test_calls_body();
    auto& muturb = sim().repo().get_field("mu_turb");

    // Get shear production field
    auto& shear_prod = sim().repo().get_field("shear_prod");

    // Check for constant value in bulk
    auto shear_bulk =
        get_val_at_kindex(shear_prod, muturb, 0, 1) / 10.0_rt / 20.0_rt;
    EXPECT_NEAR(shear_bulk, m_srate, m_tol);

    // Velocity gradients assumed in setup (init_field3)
    const amrex::Real uz_bulk = m_srate / 2.0_rt;

    // Naive answer, assumes extrapolation
    const amrex::Real s_naive = m_srate;
    // Check for different value just above wall due to BC
    auto shear_wall =
        get_val_at_kindex(shear_prod, muturb, 0, 0) / 10.0_rt / 20.0_rt;
    // Check that the result is not equal to the naive value
    EXPECT_GT(std::abs(shear_wall - s_naive), m_tol);
    // Answer that accounts for BC (Neumann)
    const amrex::Real uz_wallface_neumann =
        (1.0_rt / 3.0_rt * (uz_bulk * 1.5_rt * m_dz + 2.0_rt) +
         1.0_rt * (uz_bulk * 0.5_rt * m_dz + 2.0_rt) -
         4.0_rt / 3.0_rt * (uz_bulk * 0.5_rt * m_dz + 2.0_rt)) /
        m_dz;
    // Neumann is applied to w as well
    const amrex::Real wz_wallface_neumann = uz_wallface_neumann;
    const amrex::Real s_true = std::sqrt(
        (2.0_rt * wz_wallface_neumann * wz_wallface_neumann) +
        (2.0_rt * uz_wallface_neumann * uz_wallface_neumann));
    EXPECT_NEAR(shear_wall, s_true, m_tol);
}

TEST_F(TurbLESTestBC, test_1eqKsgs_symmetricwall)
{
    populate_parameters();
    {
        amrex::ParmParse pp("zlo");
        pp.add("type", (std::string) "symmetric_wall");
    }
    one_eq_ksgs_setup_params();
    test_calls_body();
    auto& muturb = sim().repo().get_field("mu_turb");

    // Get shear production field
    auto& shear_prod = sim().repo().get_field("shear_prod");

    // Check for constant value in bulk
    auto shear_bulk =
        get_val_at_kindex(shear_prod, muturb, 0, 1) / 10.0_rt / 20.0_rt;
    EXPECT_NEAR(shear_bulk, m_srate, m_tol);

    // Velocity gradients assumed in setup (init_field3)
    const amrex::Real uz_bulk = m_srate / 2.0_rt;

    // Naive answer, assumes extrapolation
    const amrex::Real s_naive = m_srate;
    // Check for different value just above wall due to BC
    auto shear_wall =
        get_val_at_kindex(shear_prod, muturb, 0, 0) / 10.0_rt / 20.0_rt;
    // Check that the result is not equal to the naive value
    EXPECT_GT(std::abs(shear_wall - s_naive), m_tol);
    // Answer that accounts for BC (Neumann)
    const amrex::Real uz_wallface_neumann =
        (1.0_rt / 3.0_rt * (uz_bulk * 1.5_rt * m_dz + 2.0_rt) +
         1.0_rt * (uz_bulk * 0.5_rt * m_dz + 2.0_rt) -
         4.0_rt / 3.0_rt * (uz_bulk * 0.5_rt * m_dz + 2.0_rt)) /
        m_dz;
    // Wall condition on w
    const amrex::Real wz_wallface =
        (1.0_rt / 3.0_rt * (uz_bulk * 1.5_rt * m_dz + 2.0_rt) +
         1.0_rt * (uz_bulk * 0.5_rt * m_dz + 2.0_rt) -
         4.0_rt / 3.0_rt * (0.0_rt)) /
        m_dz;
    const amrex::Real s_true = std::sqrt(
        (2.0_rt * wz_wallface * wz_wallface) +
        (2.0_rt * uz_wallface_neumann * uz_wallface_neumann));
    EXPECT_NEAR(shear_wall, s_true, m_tol);
}

} // namespace amr_wind_tests
