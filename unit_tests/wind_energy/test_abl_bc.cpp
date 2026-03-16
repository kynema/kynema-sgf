#include "abl_test_utils.H"
#include "amr-wind/utilities/trig_ops.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/incflo.H"

#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "amr-wind/equation_systems/icns/icns.H"
#include "amr-wind/equation_systems/icns/icns_ops.H"
#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/equation_systems/icns/source_terms/BodyForce.H"
#include "amr-wind/equation_systems/icns/source_terms/ABLForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/GeostrophicForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/CoriolisForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/BoussinesqBuoyancy.H"
#include "amr-wind/equation_systems/icns/source_terms/DensityBuoyancy.H"
#include "amr-wind/equation_systems/icns/source_terms/HurricaneForcing.H"
#include "amr-wind/equation_systems/icns/source_terms/RayleighDamping.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

namespace {
amrex::Real
get_val_at_kindex(amr_wind::Field& field, const int comp, const int kref)
{
    const int lev = 0;
    amrex::Real error_total = 0;

    error_total += amrex::ReduceSum(
        field(lev), 0,
        [=] AMREX_GPU_HOST_DEVICE(
            amrex::Box const& bx,
            amrex::Array4<amrex::Real const> const& f_arr) -> amrex::Real {
            amrex::Real error = 0;

            amrex::Loop(bx, [=, &error](int i, int j, int k) {
                // Check if current cell is just above lower wall
                if (k == kref) {
                    // Add field value to output
                    error += f_arr(i, j, k, comp);
                }
            });

            return error;
        });
    amrex::ParallelDescriptor::ReduceRealSum(error_total);
    return error_total;
}
void init_velocity(amr_wind::Field& fld, amrex::Real vval, int dir)
{
    const int nlevels = fld.repo().num_active_levels();

    // Initialize entire field to 0
    fld.setVal(0.0_rt);

    for (int lev = 0; lev < nlevels; ++lev) {
        fld(lev).setVal(vval, dir, 1);
    }
}
} // namespace

using ICNSFields =
    amr_wind::pde::FieldRegOp<amr_wind::pde::ICNS, amr_wind::fvm::Godunov>;

TEST_F(ABLMeshTest, abl_local_wall_model)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    constexpr amrex::Real mu = 0.01_rt;
    constexpr amrex::Real vval = 5.0_rt;
    constexpr amrex::Real dt = 0.1_rt;
    constexpr amrex::Real kappa = 0.4_rt;
    constexpr amrex::Real z0 = 0.11_rt;
    int dir = 0;
    populate_parameters();
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<int> periodic{{1, 1, 0}};
        pp.addarr("is_periodic", periodic);
    }
    {
        amrex::ParmParse pp("zlo");
        pp.add("type", (std::string) "wall_model");
    }
    {
        amrex::ParmParse pp("zhi");
        pp.add("type", (std::string) "slip_wall");
    }
    {
        amrex::ParmParse pp("incflo");
        pp.add("diffusion_type", 0);
    }
    {
        amrex::ParmParse pp("transport");
        pp.add("viscosity", mu);
    }
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", dt);
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("wall_shear_stress_type", (std::string) "local");
        pp.add("kappa", kappa);
        pp.add("surface_roughness_z0", z0);
    }
    initialize_mesh();

    // Set up solver-related routines
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().create_turbulence_model();
    sim().init_physics();

    // Specify velocity as uniform in x direction
    auto& velocity = sim().repo().get_field("velocity");
    init_velocity(velocity, vval, dir);
    // Specify density as unity
    auto& density = sim().repo().get_field("density");
    density.setVal(1.0_rt);
    // Perform post init for physics: turns on wall model
    for (auto& pp : sim().physics()) {
        pp->post_init_actions();
    }

    // Advance states to prepare for time step
    pde_mgr.advance_states();

    // Initialize icns pde
    auto& icns_eq = pde_mgr.icns();
    icns_eq.initialize();
    // Initialize viscosity
    sim().turbulence_model().update_turbulent_viscosity(
        amr_wind::FieldState::Old, DiffusionType::Crank_Nicolson);
    icns_eq.compute_mueff(amr_wind::FieldState::Old);

    // Check test setup by verifying mu
    const auto& viscosity = sim().repo().get_field("velocity_mueff");
    EXPECT_NEAR(mu, utils::field_max(viscosity), tol);
    EXPECT_NEAR(mu, utils::field_min(viscosity), tol);

    // Zero source term and convection term to focus on diffusion
    auto& src = icns_eq.fields().src_term;
    auto& adv = icns_eq.fields().conv_term;
    src.setVal(0.0_rt);
    adv.setVal(0.0_rt);

    // Calculate diffusion term
    icns_eq.compute_diffusion_term(amr_wind::FieldState::Old);
    // Setup mask_cell array to avoid errors in solve
    auto& mask_cell = sim().repo().declare_int_field("mask_cell", 1, 1);
    mask_cell.setVal(1);
    // Compute result with just diffusion term
    icns_eq.compute_predictor_rhs(DiffusionType::Explicit);

    // Get resulting velocity in first cell
    const amrex::Real vbase = get_val_at_kindex(velocity, dir, 0) / 8 / 8;

    // Calculate expected velocity after one step
    const amrex::Real dz = sim().mesh().Geom(0).CellSizeArray()[2];
    const amrex::Real zref = 0.5_rt * dz;
    const amrex::Real utau = kappa * vval / (std::log(zref / z0));
    const amrex::Real tau_wall = std::pow(utau, 2.0_rt);
    const amrex::Real vexpct = vval + (dt * (0.0_rt - tau_wall) / dz);
    EXPECT_NEAR(vexpct, vbase, tol);
}

TEST_F(ABLMeshTest, abl_donelan_wall_model)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    constexpr amrex::Real mu = 0.01_rt;
    constexpr amrex::Real vval = 35.0_rt;
    constexpr amrex::Real dt = 0.1_rt;
    constexpr amrex::Real kappa = 0.4_rt;
    constexpr amrex::Real z0 = 0.11_rt;
    int dir = 0;
    populate_parameters();
    {
        amrex::ParmParse pp("geometry");
        amrex::Vector<int> periodic{{1, 1, 0}};
        pp.addarr("is_periodic", periodic);
    }
    {
        amrex::ParmParse pp("zlo");
        pp.add("type", (std::string) "wall_model");
    }
    {
        amrex::ParmParse pp("zhi");
        pp.add("type", (std::string) "slip_wall");
    }
    {
        amrex::ParmParse pp("incflo");
        pp.add("diffusion_type", 0);
    }
    {
        amrex::ParmParse pp("transport");
        pp.add("viscosity", mu);
    }
    {
        amrex::ParmParse pp("time");
        pp.add("fixed_dt", dt);
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("wall_shear_stress_type", (std::string) "donelan");
        pp.add("kappa", kappa);
        pp.add("surface_roughness_z0", z0);
    }
    initialize_mesh();

    // Set up solver-related routines
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().create_turbulence_model();
    sim().init_physics();

    // Specify velocity as uniform in x direction
    auto& velocity = sim().repo().get_field("velocity");
    init_velocity(velocity, vval, dir);
    // Specify density as unity
    auto& density = sim().repo().get_field("density");
    density.setVal(1.0_rt);
    // Perform post init for physics: turns on wall model
    for (auto& pp : sim().physics()) {
        pp->post_init_actions();
    }

    // Advance states to prepare for time step
    pde_mgr.advance_states();

    // Initialize icns pde
    auto& icns_eq = pde_mgr.icns();
    icns_eq.initialize();
    // Initialize viscosity
    sim().turbulence_model().update_turbulent_viscosity(
        amr_wind::FieldState::Old, DiffusionType::Crank_Nicolson);
    icns_eq.compute_mueff(amr_wind::FieldState::Old);

    // Check test setup by verifying mu
    const auto& viscosity = sim().repo().get_field("velocity_mueff");
    EXPECT_NEAR(mu, utils::field_max(viscosity), tol);
    EXPECT_NEAR(mu, utils::field_min(viscosity), tol);

    // Zero source term and convection term to focus on diffusion
    auto& src = icns_eq.fields().src_term;
    auto& adv = icns_eq.fields().conv_term;
    src.setVal(0.0_rt);
    adv.setVal(0.0_rt);

    // Calculate diffusion term
    icns_eq.compute_diffusion_term(amr_wind::FieldState::Old);
    // Setup mask_cell array to avoid errors in solve
    auto& mask_cell = sim().repo().declare_int_field("mask_cell", 1, 1);
    mask_cell.setVal(1);
    // Compute result with just diffusion term
    icns_eq.compute_predictor_rhs(DiffusionType::Explicit);

    // Get resulting velocity in first cell
    const amrex::Real vbase = get_val_at_kindex(velocity, dir, 0) / 8 / 8;

    // Calculate expected velocity after one step
    const amrex::Real dz = sim().mesh().Geom(0).CellSizeArray()[2];
    const amrex::Real tau_wall = 0.0024_rt * std::pow(vval, 2.0_rt);
    const amrex::Real vexpct = vval + (dt * (0.0_rt - tau_wall) / dz);
    EXPECT_NEAR(vexpct, vbase, tol);
}

} // namespace amr_wind_tests
