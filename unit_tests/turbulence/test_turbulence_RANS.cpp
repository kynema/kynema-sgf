#include <array>
#include "gtest/gtest.h"
#include "ks_test_utils/MeshTest.H"
#include "src/turbulence/TurbulenceModel.H"
#include "ks_test_utils/test_utils.H"
#include "src/utilities/math_ops.H"

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

TEST_F(TurbRANSTest, test_1eqKrans_production_per_unit_mass)
{
    // Stable column with surface_temp_flux = 0 and every cell below
    // ABL.meso_sponge_start: KLAxell uses l = l_s(z) and Rt = 0, so
    // Cmu(Rt) = Cmu'(Rt) = Cmu and, per unit mass,
    //   shear_prod = S^2 nu_t,  buoy_prod = -nu_t N^2,  nu_t = mu_t / rho
    const amrex::Real Tref = 300.0_rt;
    const amrex::Real gravz = 9.81_rt;
    const amrex::Real rho0 = 1.2_rt;
    const amrex::Real srate = 2.0_rt;
    const amrex::Real Tgz = 0.003_rt;
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
        amrex::Vector<amrex::Real> t_hts{0.0_rt, 4000.0_rt};
        pp.addarr("temperature_heights", t_hts);
        amrex::Vector<amrex::Real> t_vals{300.0_rt, 312.0_rt};
        pp.addarr("temperature_values", t_vals);
        pp.add("surface_temp_flux", 0.0_rt);
    }
    {
        amrex::ParmParse pp("transport");
        pp.add("reference_temperature", Tref);
    }

    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();
    sim().create_turbulence_model();
    sim().turbulence_model().post_init_actions();
    auto& tmodel = sim().turbulence_model();

    auto& vel = sim().repo().get_field("velocity");
    init_strain_field(vel, srate);
    auto& dens = sim().repo().get_field("density");
    dens.setVal(rho0);
    auto& temp = sim().repo().get_field("temperature");
    init_temperature_field(temp, Tgz);
    auto& tlscale = sim().repo().get_field("turb_lscale");
    auto& tke = sim().repo().get_field("tke");

    const amrex::Real Cmu = 0.556_rt;
    const amrex::Real lambda = 30.0_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real tke_val = 0.5_rt;
    const amrex::Real nsq = gravz * Tgz / Tref;

    // Largest |computed / expected - 1| over the domain for mu_t and both
    // production terms
    auto max_errors = [&]() {
        const auto& shear_prod = sim().repo().get_field("shear_prod");
        const auto& buoy_prod = sim().repo().get_field("buoy_prod");
        const auto& muturb = sim().repo().get_field("mu_turb");
        const auto& problo = sim().mesh().Geom(0).ProbLoArray();
        const amrex::Real dz = sim().mesh().Geom(0).CellSize()[2];
        const auto& sp_arrs = shear_prod(0).const_arrays();
        const auto& bp_arrs = buoy_prod(0).const_arrays();
        const auto& mu_arrs = muturb(0).const_arrays();
        const auto r = amrex::ParReduce(
            amrex::TypeList<
                amrex::ReduceOpMax, amrex::ReduceOpMax, amrex::ReduceOpMax>{},
            amrex::TypeList<amrex::Real, amrex::Real, amrex::Real>{},
            shear_prod(0),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k)
                -> amrex::GpuTuple<amrex::Real, amrex::Real, amrex::Real> {
                const amrex::Real z = problo[2] + ((k + 0.5_rt) * dz);
                const amrex::Real ls =
                    (lambda * kappa * z) / (lambda + (kappa * z));
                const amrex::Real nut = Cmu * ls * std::sqrt(tke_val);
                const amrex::Real mu_ref = rho0 * nut;
                const amrex::Real sp_ref = srate * srate * nut;
                const amrex::Real bp_ref = -nut * nsq;
                return {
                    std::abs((mu_arrs[nbx](i, j, k) - mu_ref) / mu_ref),
                    std::abs((sp_arrs[nbx](i, j, k) - sp_ref) / sp_ref),
                    std::abs((bp_arrs[nbx](i, j, k) - bp_ref) / bp_ref)};
            });
        std::array<amrex::Real, 3> out{
            amrex::get<0>(r), amrex::get<1>(r), amrex::get<2>(r)};
        amrex::ParallelDescriptor::ReduceRealMax(out.data(), 3);
        return out;
    };

    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e3_rt;
    tke.setVal(tke_val);
    tlscale.setVal(1.0_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    auto err = max_errors();
    EXPECT_NEAR(err[0], 0.0_rt, tol);
    EXPECT_NEAR(err[1], 0.0_rt, tol);
    EXPECT_NEAR(err[2], 0.0_rt, tol);

    // Same check through the terrain kernel: flat terrain at z = 0 with no
    // blanked cells gives the same z, l_s, mu_t and production
    {
        auto& blank = sim().repo().declare_int_field("terrain_blank", 1, 1, 1);
        blank.setVal(0);
        auto& height = sim().repo().declare_field("terrain_height", 1, 1, 1);
        height.setVal(0.0_rt);
    }
    tke.setVal(tke_val);
    tlscale.setVal(1.0_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    err = max_errors();
    EXPECT_NEAR(err[0], 0.0_rt, tol);
    EXPECT_NEAR(err[1], 0.0_rt, tol);
    EXPECT_NEAR(err[2], 0.0_rt, tol);
}

} // namespace kynema_sgf_tests
