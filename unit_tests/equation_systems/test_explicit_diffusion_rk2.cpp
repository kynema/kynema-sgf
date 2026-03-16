#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/incflo_enums.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/equation_systems/PDEBase.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

namespace {

void init_scalar(amr_wind::Field& scalar)
{
    const int nlevels = scalar.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& sarrs = scalar(lev).arrays();

        amrex::ParallelFor(
            scalar(lev), scalar.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                sarrs[nbx](i, j, k) = 1.0_rt + (i * i) +
                                      (0.2_rt * (j - 1) * (j - 1)) +
                                      (0.01_rt * (k + 1) * (k + 1) * (k + 1));
            });
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace

class ExplicitDiffusionRK2Test : public MeshTest
{
protected:
    void populate_parameters() override
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
            pp.add("diffusion_type", 0);
            pp.add("density", m_rho_0);
        }
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", m_dt);
        }
    }

    void init_sim_calc_diffusion()
    {

        // High-level setup
        populate_parameters();
        initialize_mesh();

        // Initialize scalar with multiply_rho = true
        auto& pde_mgr = sim().pde_manager();
        pde_mgr.register_icns();
        sim().create_turbulence_model();
        sim().init_physics();
        auto& k_eqn = pde_mgr.register_transport_pde("TKE");
        auto& tke = k_eqn.fields().field;
        auto& tke_old = tke.state(amr_wind::FieldState::Old);
        auto& density = sim().repo().get_field("density");
        auto& density_old = density.state(amr_wind::FieldState::Old);
        auto& density_nph = density.state(amr_wind::FieldState::NPH);
        // Set up initial arrays
        density.setVal(m_rho_0);
        density_old.setVal(m_rho_0);
        density_nph.setVal(m_rho_0);
        init_scalar(tke);
        init_scalar(tke_old);
        // Zero other RHS terms for tke equation
        k_eqn.fields().conv_term.setVal(0.0_rt);
        k_eqn.fields().src_term.setVal(0.0_rt);
        // Set up mask field
        auto& mask_cell = sim().repo().declare_int_field("mask_cell", 1, 1);
        mask_cell.setVal(1);
        // Set up diffusion coefficient
        k_eqn.fields().mueff.setVal(0.1_rt);

        k_eqn.initialize();
        // Calculate diffusion term explicitly
        k_eqn.compute_diffusion_term(amr_wind::FieldState::Old);
        // Incorporate RHS, which should only be the explicit diffusion
        k_eqn.compute_predictor_rhs(DiffusionType::Explicit);
    }

    const amrex::Real m_rho_0 = 2.0_rt;
    const amrex::Real m_dt = 0.5_rt;
    const int m_nx = 8;
    const int m_ny = 8;
    const int m_nz = 16;
};

TEST_F(ExplicitDiffusionRK2Test, old_approach)
{
    init_sim_calc_diffusion();
    auto& pde_mgr = sim().pde_manager();
    auto& k_eqn = pde_mgr("TKE-Godunov");
    auto& tke = k_eqn.fields().field;
    auto& tke_old = tke.state(amr_wind::FieldState::Old);

    // Subtract diffusion term manually, mimicking the assumption in the code
    // If correct, should go back to old value
    const auto dt = sim().time().delta_t();
    auto& diff_old = k_eqn.fields().diff_term.state(amr_wind::FieldState::New);
    amr_wind::field_ops::saxpy(tke, -dt, diff_old, 0, 0, 1, 0);

    // Subtract original tke from result for comparison
    amr_wind::field_ops::saxpy(tke, -1.0_rt, tke_old, 0, 0, 1, 0);

    auto min_diff = utils::field_min(tke);
    auto max_diff = utils::field_max(tke);
    const auto abs_diff_tke = amrex::max<amrex::Real>(-min_diff, max_diff);
    // We know this approach is wrong, so expect error to be greater
    EXPECT_GT(abs_diff_tke, std::numeric_limits<float>::epsilon() * 1.0e1_rt);

    // As a sanity check for the test, confirm that diff_term is not 0
    min_diff = utils::field_min(diff_old);
    max_diff = utils::field_max(diff_old);
    const auto abs_diff_term = amrex::max<amrex::Real>(-min_diff, max_diff);
    EXPECT_GT(abs_diff_term, std::numeric_limits<float>::epsilon() * 1.0e1_rt);
}

TEST_F(ExplicitDiffusionRK2Test, correct_approach)
{
    init_sim_calc_diffusion();
    auto& pde_mgr = sim().pde_manager();
    auto& k_eqn = pde_mgr("TKE-Godunov");
    auto& tke = k_eqn.fields().field;
    auto& tke_old = tke.state(amr_wind::FieldState::Old);

    // Change term to -2 times original value so it will cancel RHS contribution
    const auto dt = sim().time().delta_t();
    auto& diff_new = k_eqn.fields().diff_term.state(amr_wind::FieldState::New);
    amr_wind::field_ops::saxpy(diff_new, -3.0_rt, diff_new, 0, 0, 1, 0);

    // Apply operation through PDE function
    k_eqn.improve_explicit_diffusion(dt);

    // Subtract original tke from result for comparison
    amr_wind::field_ops::saxpy(tke, -1.0_rt, tke_old, 0, 0, 1, 0);

    // Check that difference is less than tolerance
    auto min_diff = utils::field_min(tke);
    auto max_diff = utils::field_max(tke);
    const auto abs_diff_tke = amrex::max<amrex::Real>(-min_diff, max_diff);
    EXPECT_LT(
        abs_diff_tke, std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);

    // As a sanity check for the test, confirm that diff_term is not 0
    min_diff = utils::field_min(diff_new);
    max_diff = utils::field_max(diff_new);
    const auto abs_diff_term = amrex::max<amrex::Real>(-min_diff, max_diff);
    EXPECT_GT(
        abs_diff_term, std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
}

} // namespace amr_wind_tests
