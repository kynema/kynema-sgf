#include "gtest/gtest.h"
#include "ks_test_utils/MeshTest.H"
#include "src/turbulence/TurbulenceModel.H"
#include "ks_test_utils/test_utils.H"
#include "src/utilities/math_ops.H"
#include "src/equation_systems/tke/source_terms/KransAxell.H"
#include <fstream>

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

void init_strain_field(kynema_sgf::Field& fld, amrex::Real srate)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();
    amrex::Real offset =
        (fld.field_location() == kynema_sgf::FieldLoc::CELL) ? 0.5_rt : 0.0_rt;
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = fld(lev).arrays();

        amrex::ParallelFor(
            fld(lev), fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real x = problo[0] + ((i + offset) * dx[0]);
                const amrex::Real y = problo[1] + ((j + offset) * dx[1]);
                const amrex::Real z = problo[2] + ((k + offset) * dx[2]);

                farrs[nbx](i, j, k, 0) = x / std::sqrt(6.0_rt) * srate;
                farrs[nbx](i, j, k, 1) = y / std::sqrt(6.0_rt) * srate;
                farrs[nbx](i, j, k, 2) = z / std::sqrt(6.0_rt) * srate;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void init_temperature_field(kynema_sgf::Field& fld, amrex::Real tgrad)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();

    amrex::Real offset =
        (fld.field_location() == kynema_sgf::FieldLoc::CELL) ? 0.5_rt : 0.0_rt;

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

// value below ztop, zero above
void set_below(
    kynema_sgf::Field& fld, const amrex::Real value, const amrex::Real ztop)
{
    const auto& mesh = fld.repo().mesh();
    for (int lev = 0; lev < fld.repo().num_active_levels(); ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = fld(lev).arrays();
        amrex::ParallelFor(
            fld(lev), fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
                farrs[nbx](i, j, k) = (z < ztop) ? value : 0.0_rt;
            });
    }
    amrex::Gpu::streamSynchronize();
}

// Neutral KLAxell inputs of test_1eqKrans_setup_calc, plus the 1-D RANS
// profile needed by the KransAxell source
void add_klaxell_inputs(
    const amrex::Real Tref, const amrex::Real gravz, const amrex::Real rho0)
{
    const std::string rans_file = "rans_1d_canopy.info";
    {
        std::ofstream os(rans_file);
        os << "0 8 0 0 0.1\n";
        os << "500 8 0 0 0.1\n";
        os << "4000 8 0 0 0.1\n";
    }
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "KLAxell");
    }
    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<std::string> physics{"ABL"};
        pp.addarr("physics", physics);
        pp.add("density", rho0);
        amrex::Vector<amrex::Real> vvec{8.0_rt, 0.0_rt, 0.0_rt};
        pp.addarr("velocity", vvec);
        amrex::Vector<amrex::Real> gvec{0.0_rt, 0.0_rt, -gravz};
        pp.addarr("gravity", gvec);
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("surface_temp_rate", 0.0_rt);
        pp.add("initial_wind_profile", true);
        amrex::Vector<amrex::Real> t_hts{0.0_rt, 100.0_rt, 4000.0_rt};
        pp.addarr("temperature_heights", t_hts);
        pp.addarr("wind_heights", t_hts);
        amrex::Vector<amrex::Real> t_vals{Tref, Tref, Tref};
        pp.addarr("temperature_values", t_vals);
        amrex::Vector<amrex::Real> u_vals{8.0_rt, 8.0_rt, 8.0_rt};
        pp.addarr("u_values", u_vals);
        amrex::Vector<amrex::Real> v_vals{0.0_rt, 0.0_rt, 0.0_rt};
        pp.addarr("v_values", v_vals);
        amrex::Vector<amrex::Real> tke_vals{0.1_rt, 0.1_rt, 0.1_rt};
        pp.addarr("tke_values", tke_vals);
        pp.add("surface_temp_flux", 0.0_rt);
        pp.add("rans_1dprofile_file", rans_file);
    }
    {
        amrex::ParmParse pp("transport");
        pp.add("reference_temperature", Tref);
    }
}

} // namespace

class TurbRANSTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{10, 10, 64}};
            pp.addarr("n_cell", ncell);
            pp.add("blocking_factor", 2);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{
                {1024.0_rt, 1024.0_rt, 1024.0_rt}};
            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
    }
};

TEST_F(TurbRANSTest, test_1eqKrans_setup_calc)
{
    // Parser inputs for turbulence model
    const amrex::Real Tref = 265.0_rt;
    const amrex::Real gravz = 10.0_rt;
    const amrex::Real rho0 = 1.2_rt;
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "KLAxell");
    }
    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<std::string> physics{"ABL"};
        pp.addarr("physics", physics);
        pp.add("density", rho0);
        amrex::Vector<amrex::Real> vvec{8.0_rt, 0.0_rt, 0.0_rt};
        pp.addarr("velocity", vvec);
        amrex::Vector<amrex::Real> gvec{0.0_rt, 0.0_rt, -gravz};
        pp.addarr("gravity", gvec);
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("surface_temp_rate", 0.0_rt);
        pp.add("initial_wind_profile", true);
        amrex::Vector<amrex::Real> t_hts{0.0_rt, 100.0_rt, 4000.0_rt};
        pp.addarr("temperature_heights", t_hts);
        pp.addarr("wind_heights", t_hts);
        amrex::Vector<amrex::Real> t_vals{265.0_rt, 265.0_rt, 265.0_rt};
        pp.addarr("temperature_values", t_vals);
        amrex::Vector<amrex::Real> u_vals{8.0_rt, 8.0_rt, 8.0_rt};
        pp.addarr("u_values", u_vals);
        amrex::Vector<amrex::Real> v_vals{0.0_rt, 0.0_rt, 0.0_rt};
        pp.addarr("v_values", v_vals);
        amrex::Vector<amrex::Real> tke_vals{0.1_rt, 0.1_rt, 0.1_rt};
        pp.addarr("tke_values", tke_vals);
        pp.add("surface_temp_flux", 0.0_rt);
    }
    // Transport
    {
        amrex::ParmParse pp("transport");
        pp.add("reference_temperature", Tref);
    }

    // Initialize necessary parts of solver
    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();

    // Create turbulence model
    sim().create_turbulence_model();
    sim().turbulence_model().post_init_actions();
    // Get turbulence model
    auto& tmodel = sim().turbulence_model();

    // Get coefficients
    auto model_dict = tmodel.model_coeffs();

    // Constants for fields
    const amrex::Real srate = 20.0_rt;
    const amrex::Real Tgz = 0.0_rt;
    const amrex::Real lambda = 30.0_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real x3 = 1016.0_rt;
    const amrex::Real lscale_s =
        (lambda * kappa * x3) / (lambda + (kappa * x3));
    const amrex::Real tlscale_val = lscale_s;
    const amrex::Real tke_val = 0.1_rt;
    // Set up velocity field with constant strainrate
    auto& vel = sim().repo().get_field("velocity");
    init_strain_field(vel, srate);
    // Set up uniform unity density field
    auto& dens = sim().repo().get_field("density");
    dens.setVal(rho0);
    // Set up temperature field with constant gradient in z
    auto& temp = sim().repo().get_field("temperature");
    init_temperature_field(temp, Tgz);
    // Give values to tlscale and tke arrays
    auto& tlscale = sim().repo().get_field("turb_lscale");
    tlscale.setVal(tlscale_val);
    auto& tke = sim().repo().get_field("tke");
    tke.setVal(tke_val);

    // Update turbulent viscosity directly
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const auto& muturb = sim().repo().get_field("mu_turb");

    // Check values of turbulent viscosity
    const auto max_val = utils::field_max(muturb);
    const amrex::Real Cmu = 0.556_rt;
    const amrex::Real epsilon = kynema_sgf::utils::powi(Cmu, 3) *
                                std::pow(tke_val, 1.5_rt) /
                                (tlscale_val + 1.0e-3_rt);
    const amrex::Real stratification = 0.0_rt;
    const amrex::Real Rt =
        kynema_sgf::utils::powi(tke_val / epsilon, 2) * stratification;
    const amrex::Real Cmu_Rt = (0.556_rt + (0.108_rt * Rt)) /
                               (1.0_rt + (0.308_rt * Rt) +
                                (0.00837_rt * kynema_sgf::utils::powi(Rt, 2)));
    const amrex::Real tol = 0.12_rt;
    const amrex::Real nut_max =
        rho0 * Cmu_Rt * tlscale_val * std::sqrt(tke_val);
    EXPECT_NEAR(max_val, nut_max, tol);
}

// ForestDrag.canopy_tke adds fd (beta_p |U|^3 - beta_d |U| k) to the
// KransAxell source in the canopy only; off, a forest changes nothing
TEST_F(TurbRANSTest, test_1eqKrans_canopy_tke_source)
{
    const amrex::Real rho0 = 1.2_rt;
    add_klaxell_inputs(265.0_rt, 10.0_rt, rho0);
    populate_parameters();
    initialize_mesh();
    sim().pde_manager().register_icns();
    sim().init_physics();
    sim().create_transport_model();
    sim().create_turbulence_model();
    sim().turbulence_model().post_init_actions();
    auto& tmodel = sim().turbulence_model();

    auto& repo = sim().repo();
    auto& vel = repo.get_field("velocity");
    init_strain_field(vel, 0.01_rt);
    repo.get_field("density").setVal(rho0);
    init_temperature_field(repo.get_field("temperature"), 0.0_rt);
    repo.get_field("turb_lscale").setVal(10.0_rt);
    const amrex::Real tke_val = 0.1_rt;
    repo.get_field("tke").setVal(tke_val);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const amrex::Real dt = 0.5_rt;
    sim().time().delta_t() = dt;

    auto& src_no_forest = repo.declare_field("src_no_forest", 1, 0, 1);
    auto& src_off = repo.declare_field("src_off", 1, 0, 1);
    auto& src_on = repo.declare_field("src_on", 1, 0, 1);
    src_no_forest.setVal(0.0_rt);
    src_off.setVal(0.0_rt);
    src_on.setVal(0.0_rt);

    // Default inputs, before the forest exists
    kynema_sgf::pde::tke::KransAxell krans_off(sim());
    krans_off(0, kynema_sgf::FieldState::New, src_no_forest(0));

    // Canopy below 100 m
    const amrex::Real fd = 0.02_rt;
    auto& forest = repo.declare_field("forest_drag", 1, 1, 1);
    set_below(forest, fd, 100.0_rt);
    krans_off(0, kynema_sgf::FieldState::New, src_off(0));

    {
        amrex::ParmParse pp("ForestDrag");
        pp.add("canopy_tke", true);
    }
    kynema_sgf::pde::tke::KransAxell krans_on(sim());
    krans_on(0, kynema_sgf::FieldState::New, src_on(0));

    // Off: bit-identical to no forest
    amrex::MultiFab::Subtract(src_off(0), src_no_forest(0), 0, 0, 1, 0);
    EXPECT_EQ(src_off(0).norm0(), 0.0_rt);

    // On: the canopy term in the canopy (z = 88 m), nothing above (z = 104 m)
    amrex::MultiFab::Subtract(src_on(0), src_no_forest(0), 0, 0, 1, 0);
    const int i = 5;
    const int j = 5;
    const int kc = 5;
    const amrex::Real ux = utils::field_probe(vel, 0, i, j, kc, 0);
    const amrex::Real uy = utils::field_probe(vel, 0, i, j, kc, 1);
    const amrex::Real uz = utils::field_probe(vel, 0, i, j, kc, 2);
    const amrex::Real ws = std::sqrt((ux * ux) + (uy * uy) + (uz * uz));
    const amrex::Real expected =
        (fd * ws * ws * ws) +
        (tke_val * std::expm1(-4.0_rt * fd * ws * dt) / dt);
    EXPECT_NEAR(
        utils::field_probe(src_on, 0, i, j, kc), expected,
        1.0e-6_rt * std::abs(expected));
    EXPECT_EQ(utils::field_probe(src_on, 0, i, j, kc + 1), 0.0_rt);
    EXPECT_EQ(utils::field_probe(src_on, 0, i, j, 40), 0.0_rt);
}

// KLAxell with the forest present and ForestDrag.canopy_tke off: the eddy
// viscosity is bit-identical to the run without a forest
TEST_F(TurbRANSTest, test_1eqKrans_canopy_length_off)
{
    const amrex::Real rho0 = 1.2_rt;
    add_klaxell_inputs(265.0_rt, 10.0_rt, rho0);
    populate_parameters();
    initialize_mesh();
    sim().pde_manager().register_icns();
    sim().init_physics();
    sim().create_transport_model();
    sim().create_turbulence_model();
    sim().turbulence_model().post_init_actions();
    auto& tmodel = sim().turbulence_model();

    auto& repo = sim().repo();
    init_strain_field(repo.get_field("velocity"), 0.01_rt);
    repo.get_field("density").setVal(rho0);
    init_temperature_field(repo.get_field("temperature"), 0.0_rt);
    repo.get_field("turb_lscale").setVal(10.0_rt);
    repo.get_field("tke").setVal(0.1_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    auto& mu_no_forest = repo.declare_field("mu_no_forest", 1, 0, 1);
    amrex::MultiFab::Copy(
        mu_no_forest(0), repo.get_field("mu_turb")(0), 0, 0, 1, 0);

    set_below(repo.declare_field("forest_drag", 1, 1, 1), 0.02_rt, 100.0_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    amrex::MultiFab::Subtract(
        mu_no_forest(0), repo.get_field("mu_turb")(0), 0, 0, 1, 0);
    EXPECT_EQ(mu_no_forest(0).norm0(), 0.0_rt);
}

// ForestDrag.canopy_tke on: l = min(l_s, alpha / fd) in the canopy only
TEST_F(TurbRANSTest, test_1eqKrans_canopy_length_on)
{
    const amrex::Real rho0 = 1.2_rt;
    add_klaxell_inputs(265.0_rt, 10.0_rt, rho0);
    {
        amrex::ParmParse pp("ForestDrag");
        pp.add("canopy_tke", true);
    }
    populate_parameters();
    initialize_mesh();
    sim().pde_manager().register_icns();
    sim().init_physics();
    sim().create_transport_model();
    sim().create_turbulence_model();
    sim().turbulence_model().post_init_actions();
    auto& tmodel = sim().turbulence_model();

    auto& repo = sim().repo();
    init_strain_field(repo.get_field("velocity"), 0.01_rt);
    repo.get_field("density").setVal(rho0);
    init_temperature_field(repo.get_field("temperature"), 0.0_rt);
    repo.get_field("turb_lscale").setVal(10.0_rt);
    const amrex::Real tke_val = 0.1_rt;
    repo.get_field("tke").setVal(tke_val);

    // Zero forest drag: no limit anywhere
    auto& forest = repo.declare_field("forest_drag", 1, 1, 1);
    forest.setVal(0.0_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    auto& mu_free = repo.declare_field("mu_free", 1, 0, 1);
    amrex::MultiFab::Copy(mu_free(0), repo.get_field("mu_turb")(0), 0, 0, 1, 0);

    // Canopy below 100 m: alpha / fd = 2 m, below l_s at every height
    const amrex::Real fd = 0.02_rt;
    set_below(forest, fd, 100.0_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);

    const auto& mu = repo.get_field("mu_turb");
    const auto& tlscale = repo.get_field("turb_lscale");
    const amrex::Real Cmu = 0.556_rt;
    const amrex::Real lcap = 0.04_rt / fd;
    const amrex::Real mu_cap = rho0 * Cmu * lcap * std::sqrt(tke_val);
    const amrex::Real tol = 10.0_rt * kynema_sgf::constants::TIGHT_TOL;
    for (int k = 0; k < 6; ++k) {
        EXPECT_NEAR(utils::field_probe(tlscale, 0, 5, 5, k), lcap, tol);
        EXPECT_NEAR(utils::field_probe(mu, 0, 5, 5, k), mu_cap, tol);
    }
    // Above the canopy (z = 104 m and up) the eddy viscosity is unchanged
    for (int k = 6; k < 64; k += 7) {
        EXPECT_EQ(
            utils::field_probe(mu, 0, 5, 5, k),
            utils::field_probe(mu_free, 0, 5, 5, k));
    }
}

} // namespace kynema_sgf_tests
