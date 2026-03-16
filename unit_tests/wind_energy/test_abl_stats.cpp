#include "abl_test_utils.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/equation_systems/tke/TKE.H"
#include "amr-wind/wind_energy/ABLStats.H"
#include "amr-wind/incflo.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

namespace {
void init_field1(amr_wind::Field& fld)
{
    const int nlevels = fld.repo().num_active_levels();

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& farrs = fld(lev).arrays();
        amrex::ParallelFor(
            fld(lev), amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                farrs[nbx](i, j, k) =
                    std::sin(0.01_rt * i) +
                    std::pow(static_cast<amrex::Real>(k), 0.2_rt) +
                    std::cos(0.01_rt * j);
            });
    }
    amrex::Gpu::streamSynchronize();
}

amrex::Real test_new_tke(
    amr_wind::Field& tke,
    amr_wind::Field& tkeold,
    amr_wind::Field& conv_term,
    amr_wind::Field& buoy_prod,
    amr_wind::Field& shear_prod,
    amr_wind::Field& dissipation,
    amr_wind::ScratchField& diffusion,
    const amrex::Real dt)
{
    amrex::Real error_total = 0;

    for (int lev = 0; lev < tke.repo().num_active_levels(); ++lev) {

        // Form tke estimate by adding to the old tke field
        const auto& tke_old_arrs = tkeold(lev).arrays();
        const auto& buoy_prod_arrs = buoy_prod(lev).const_arrays();
        const auto& shear_prod_arrs = shear_prod(lev).const_arrays();
        const auto& dissipation_arrs = dissipation(lev).const_arrays();
        const auto& diffusion_arrs = diffusion(lev).const_arrays();
        const auto& conv_arrs = conv_term(lev).const_arrays();
        amrex::ParallelFor(
            tkeold(lev), amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                tke_old_arrs[nbx](i, j, k) +=
                    dt *
                    (conv_arrs[nbx](i, j, k) + shear_prod_arrs[nbx](i, j, k) +
                     buoy_prod_arrs[nbx](i, j, k) -
                     dissipation_arrs[nbx](i, j, k) +
                     diffusion_arrs[nbx](i, j, k));
            });
        amrex::Gpu::streamSynchronize();

        // Difference between tke estimate and tke calculated by the code
        error_total += amrex::ReduceSum(
            tke(lev), tkeold(lev), 0,
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Box const& bx,
                amrex::Array4<amrex::Real const> const& tke_arr,
                amrex::Array4<amrex::Real const> const& tke_est)
                -> amrex::Real {
                amrex::Real error = 0;

                amrex::Loop(bx, [=, &error](int i, int j, int k) {
                    error += std::abs(tke_arr(i, j, k) - tke_est(i, j, k));
                });

                return error;
            });
    }
    return error_total;
}

void remove_nans(amr_wind::Field& field)
{
    for (int lev = 0; lev < field.repo().num_active_levels(); ++lev) {
        const auto& farrs = field(lev).arrays();
        amrex::ParallelFor(
            field(lev), amrex::IntVect(0), field(lev).nComp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                farrs[nbx](i, j, k, n) = std::isnan(farrs[nbx](i, j, k, n))
                                             ? 0.0_rt
                                             : farrs[nbx](i, j, k, n);
            });
    }
    amrex::Gpu::streamSynchronize();
}
} // namespace

TEST_F(ABLMeshTest, stats_tke_diffusion)
{
    // This test checks the implementation of the calc_diffusion routine to see
    // if it does what it is intended to do
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    constexpr amrex::Real val_shear = 0.5_rt;
    constexpr amrex::Real val_buoy = 0.4_rt;
    constexpr amrex::Real val_dissip = 0.3_rt;
    constexpr amrex::Real val_conv = 0.2_rt;
    constexpr amrex::Real val_tke = 1.6_rt;
    constexpr amrex::Real val_tkeold = 1.0_rt;
    constexpr amrex::Real dt = 0.1_rt;
    constexpr amrex::Real expected_diff = ((val_tke - val_tkeold) / dt) -
                                          val_conv - val_shear - val_buoy +
                                          val_dissip;

    populate_parameters();
    initialize_mesh();

    // Register fields and eqs for the sake of a functional test
    sim().repo().declare_field("temperature", 1, 3);
    sim().pde_manager().register_icns();

    // Register fields used in turbulence model OneEqKsgs
    auto& shear = sim().repo().declare_field("shear_prod", 1, 1);
    auto& buoy = sim().repo().declare_field("buoy_prod", 1, 1);
    auto& dissip = sim().repo().declare_field("dissipation", 1, 1);
    auto diff = sim().repo().create_scratch_field("diffusion", 1, 1);
    // Register transport pde for tke
    auto& tke_eqn = sim().pde_manager().register_transport_pde(
        amr_wind::pde::TKE::pde_name());
    auto& tke = tke_eqn.fields().field;

    // Populate values of fields
    shear.setVal(val_shear);
    buoy.setVal(val_buoy);
    dissip.setVal(val_dissip);
    tke_eqn.fields().conv_term.setVal(val_conv);
    tke.setVal(val_tke);
    tke.state(amr_wind::FieldState::Old).setVal(val_tkeold);

    // Initialize ABL Stats
    amr_wind::ABLWallFunction wall_func(sim());
    amr_wind::ABLStats stats(sim(), wall_func, 2);

    // Calculate diffusion term
    stats.calc_tke_diffusion(*diff, buoy, shear, dissip, dt);

    // Check answer
    const auto min_val = utils::field_min(*diff, 0);
    const auto max_val = utils::field_max(*diff, 0);
    EXPECT_NEAR(expected_diff, min_val, tol);
    EXPECT_NEAR(expected_diff, max_val, tol);
}

TEST_F(ABLMeshTest, stats_energy_budget)
{
    // This test checks the assumptions behind the calc_diffusion routine
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    constexpr amrex::Real dt = 0.1_rt;
    constexpr amrex::Real val_shear = 0.5_rt;
    constexpr amrex::Real val_buoy = 0.4_rt;
    constexpr amrex::Real val_tlscale = 0.3_rt;

    // Set up turbulence model and other input arguments
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "OneEqKsgsM84");
    }
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", dt);
    }
    {
        amrex::ParmParse pp("TKE");
        pp.add("source_terms", (std::string) "KsgsM84Src");
    }
    {
        amrex::ParmParse pp("transport");
        pp.add("viscosity", 1.0e-5_rt);
    }

    // incflo.diffusion_type = 1
    populate_parameters();
    initialize_mesh();

    // Register icns
    auto& icns = sim().pde_manager().register_icns();
    icns.initialize();

    // Initialize fields
    sim().init_physics();

    // Initialize turbulence model
    sim().create_turbulence_model();

    // Initialize ABL velocity
    const int lev = 0;
    for (auto& pp : sim().physics()) {
        pp->pre_init_actions();
        pp->initialize_fields(lev, sim().mesh().Geom(lev));
    }

    // Register transport pde for tke
    auto& tke_eqn = sim().pde_manager().register_transport_pde(
        amr_wind::pde::TKE::pde_name());
    tke_eqn.initialize();
    auto& tke = tke_eqn.fields().field;

    // Initialize ABL Stats
    amr_wind::ABLWallFunction wall_func(sim());
    amr_wind::ABLStats stats(sim(), wall_func, 2);

    // Set initial tke value and advance states
    init_field1(tke);
    sim().pde_manager().advance_states();

    // Initialize fields for tke source term
    auto& shear = sim().repo().get_field("shear_prod");
    auto& buoy = sim().repo().get_field("buoy_prod");
    auto& tlscale = sim().repo().get_field("turb_lscale");
    shear.setVal(val_shear);
    buoy.setVal(val_buoy);
    tlscale.setVal(val_tlscale);

    // Set up new and NPH density for the sake of src term
    auto& density_nph =
        sim().repo().get_field("density").state(amr_wind::FieldState::NPH);
    density_nph.setVal(1.0_rt);

    // Setup mask_cell array to avoid errors in solve
    auto& mask_cell = sim().repo().declare_int_field("mask_cell", 1, 1);
    mask_cell.setVal(1);

    // Populate advection velocities
    icns.pre_advection_actions(amr_wind::FieldState::Old);

    // Step forward in time for tke equation
    sim().turbulence_model().update_turbulent_viscosity(
        amr_wind::FieldState::Old, DiffusionType::Crank_Nicolson);
    tke_eqn.compute_advection_term(amr_wind::FieldState::Old);
    // Remove NaNs (not sure why they're there, but need to be removed)
    remove_nans(tke_eqn.fields().conv_term);
    tke_eqn.compute_mueff(amr_wind::FieldState::Old);
    tke_eqn.compute_source_term(amr_wind::FieldState::NPH);
    tke_eqn.compute_diffusion_term(amr_wind::FieldState::New);
    tke_eqn.compute_predictor_rhs(DiffusionType::Crank_Nicolson);
    tke_eqn.solve(0.5_rt * dt);
    tke_eqn.post_solve_actions();

    // Calculate diffusion term
    auto& dissip = sim().repo().get_field("dissipation");
    auto diff = sim().repo().create_scratch_field("diffusion", 1, 1);
    stats.calc_tke_diffusion(*diff, buoy, shear, dissip, dt);

    // Check assumptions in diffusion term: sum of terms gets result
    const amrex::Real err_total = test_new_tke(
        tke, tke.state(amr_wind::FieldState::Old), tke_eqn.fields().conv_term,
        buoy, shear, dissip, *(diff), dt);
    EXPECT_NEAR(err_total, 0.0_rt, tol);
}

} // namespace amr_wind_tests
