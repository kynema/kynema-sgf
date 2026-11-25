#include <array>
#include "gtest/gtest.h"
#include "ks_test_utils/MeshTest.H"
#include "src/core/field_ops.H"
#include "src/turbulence/TurbulenceModel.H"
#include "ks_test_utils/test_utils.H"
#include "src/utilities/constants.H"
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

void init_linear_z(kynema_sgf::Field& fld, amrex::Real base, amrex::Real grad)
{
    const auto& mesh = fld.repo().mesh();
    for (int lev = 0; lev < fld.repo().num_active_levels(); ++lev) {
        const auto& geom = mesh.Geom(lev);
        const amrex::Real zlo = geom.ProbLoArray()[2];
        const amrex::Real dz = geom.CellSizeArray()[2];
        const auto& arrs = fld(lev).arrays();
        amrex::ParallelFor(
            fld(lev), fld.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                arrs[nbx](i, j, k) = base + (grad * (zlo + (k + 0.5_rt) * dz));
            });
    }
    amrex::Gpu::streamSynchronize();
}

std::array<amrex::Real, 5> kw_source_errors(
    kynema_sgf::Field& shear_prod,
    amrex::Real srate,
    bool check_buoyancy,
    bool check_cross_diffusion)
{
    const auto& repo = shear_prod.repo();
    const auto& rho = repo.get_field("density")(0).const_arrays();
    const auto& tke = repo.get_field("tke")(0).const_arrays();
    const auto& sdr = repo.get_field("sdr")(0).const_arrays();
    const auto& mu = repo.get_field("mu_turb")(0).const_arrays();
    const auto& f1 = repo.get_field("f1")(0).const_arrays();
    const auto& shear = shear_prod(0).const_arrays();
    const auto& buoy = repo.get_field("buoyancy_term")(0).const_arrays();
    const auto& cross = repo.get_field("sdr_dissipation")(0).const_arrays();
    const auto result = amrex::ParReduce(
        amrex::TypeList<
            amrex::ReduceOpMax, amrex::ReduceOpMax, amrex::ReduceOpMax,
            amrex::ReduceOpMin, amrex::ReduceOpMax>{},
        amrex::TypeList<
            amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real>{},
        shear_prod(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k)
            -> amrex::GpuTuple<
                amrex::Real, amrex::Real, amrex::Real, amrex::Real,
                amrex::Real> {
            const amrex::Real density = rho[nbx](i, j, k);
            const amrex::Real uncapped =
                mu[nbx](i, j, k) / density * srate * srate;
            const amrex::Real limit =
                0.9_rt * tke[nbx](i, j, k) * sdr[nbx](i, j, k);
            const amrex::Real expected_buoy = 9.81_rt * 0.0002_rt *
                                              mu[nbx](i, j, k) /
                                              (density * density * 0.85_rt);
            const amrex::Real expected_cross =
                (1.0_rt - f1[nbx](i, j, k)) * 2.0_rt * 0.856_rt * 0.0001_rt *
                0.0002_rt /
                (sdr[nbx](i, j, k) +
                 std::numeric_limits<amrex::Real>::epsilon() * 1.0e1_rt);
            return {
                std::abs(
                    shear[nbx](i, j, k) -
                    amrex::min<amrex::Real>(uncapped, limit)),
                check_buoyancy ? std::abs(buoy[nbx](i, j, k) - expected_buoy)
                               : 0.0_rt,
                check_cross_diffusion
                    ? std::abs(cross[nbx](i, j, k) - expected_cross)
                    : 0.0_rt,
                uncapped - limit, uncapped - limit};
        });
    std::array<amrex::Real, 5> out{
        amrex::get<0>(result), amrex::get<1>(result), amrex::get<2>(result),
        amrex::get<3>(result), amrex::get<4>(result)};
    amrex::ParallelDescriptor::ReduceRealMax(out.data(), 3);
    amrex::ParallelDescriptor::ReduceRealMin(out.data() + 3, 1);
    amrex::ParallelDescriptor::ReduceRealMax(out.data() + 4, 1);
    return out;
}

// Largest |computed / expected - 1| for mu_t and both production terms
std::array<amrex::Real, 3> max_production_errors(
    kynema_sgf::Field& shear_prod,
    amrex::Real Cmu,
    amrex::Real lambda,
    amrex::Real kappa,
    amrex::Real tke_val,
    amrex::Real rho0,
    amrex::Real srate,
    amrex::Real nsq)
{
    const auto& buoy_prod = shear_prod.repo().get_field("buoy_prod");
    const auto& muturb = shear_prod.repo().get_field("mu_turb");
    const auto& geom = shear_prod.repo().mesh().Geom(0);
    const auto& problo = geom.ProbLoArray();
    const amrex::Real dz = geom.CellSize()[2];
    const auto& sp_arrs = shear_prod(0).const_arrays();
    const auto& bp_arrs = buoy_prod(0).const_arrays();
    const auto& mu_arrs = muturb(0).const_arrays();
    const auto r = amrex::ParReduce(
        amrex::TypeList<
            amrex::ReduceOpMax, amrex::ReduceOpMax, amrex::ReduceOpMax>{},
        amrex::TypeList<amrex::Real, amrex::Real, amrex::Real>{}, shear_prod(0),
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
}

} // namespace

class TurbRANSTest : public MeshTest
{
protected:
    void check_kw_density_sources(const std::string& model);

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

void TurbRANSTest::check_kw_density_sources(const std::string& model)
{
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", model);
    }
    {
        amrex::ParmParse pp("incflo");
        pp.add("constant_density", false);
    }
    {
        amrex::ParmParse pp("transport");
        pp.add("viscosity", 0.0_rt);
    }
    {
        amrex::ParmParse pp(model + "_coeffs");
        pp.add("tke_amb", 0.1_rt);
        pp.add("sdr_amb", 0.2_rt);
    }
    {
        amrex::ParmParse pp("SDR");
        amrex::Vector<std::string> sources{"SDRSrc"};
        pp.addarr("source_terms", sources);
    }

    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().init_physics();
    sim().create_turbulence_model();
    auto& tmodel = sim().turbulence_model();
    tmodel.post_init_actions();
    auto& tke_eqn = pde_mgr("TKE-" + pde_mgr.scheme());
    auto& sdr_eqn = pde_mgr("SDR-" + pde_mgr.scheme());
    tke_eqn.initialize();
    sdr_eqn.initialize();
    time().delta_t() = 0.1_rt;

    auto& repo = sim().repo();
    auto& rho = repo.get_field("density");
    auto& rho_nph = rho.state(kynema_sgf::FieldState::NPH);
    auto& vel = repo.get_field("velocity");
    init_linear_z(repo.get_field("tke"), 0.5_rt, 0.0001_rt);
    init_linear_z(repo.get_field("sdr"), 0.5_rt, 0.0002_rt);
    repo.get_field("wall_dist").setVal(20.0_rt);

    const std::array<std::pair<std::string, amrex::Real>, 9> fields{{
        {"shear_prod", 1.0_rt},
        {"dissipation", 1.0_rt},
        {"buoyancy_term", 1.0_rt},
        {"omega_src", 1.0_rt},
        {"sdr_dissipation", 1.0_rt},
        {"tke_lhs_src_term", 2.5_rt},
        {"sdr_lhs_src_term", 2.5_rt},
        {"tke_src_term", 2.5_rt},
        {"sdr_src_term", 2.5_rt},
    }};
    const std::array<DiffusionType, 3> diffusion_types{
        DiffusionType::Explicit, DiffusionType::Crank_Nicolson,
        DiffusionType::Implicit};

    for (const auto diff_type : diffusion_types) {
        for (const amrex::Real srate : {0.2_rt, 2.0_rt}) {
            SCOPED_TRACE(model);
            SCOPED_TRACE(static_cast<int>(diff_type));
            SCOPED_TRACE(srate);
            init_strain_field(vel, srate);
            init_linear_z(rho, 1.0_rt, 0.0002_rt);
            init_linear_z(rho_nph, 1.0_rt, 0.0002_rt);
            tmodel.update_turbulent_viscosity(
                kynema_sgf::FieldState::New, diff_type);
            tke_eqn.compute_source_term(kynema_sgf::FieldState::NPH);
            sdr_eqn.compute_source_term(kynema_sgf::FieldState::NPH);

            const auto errors = kw_source_errors(
                repo.get_field("shear_prod"), srate, model == "KOmegaSST",
                diff_type == DiffusionType::Crank_Nicolson);
            const amrex::Real reference_tol = kynema_sgf::constants::TIGHT_TOL;
            EXPECT_NEAR(
                errors[0], 0.0_rt,
                reference_tol * utils::field_max(repo.get_field("shear_prod")));
            EXPECT_NEAR(
                errors[1], 0.0_rt,
                reference_tol *
                    std::abs(
                        utils::field_max(repo.get_field("buoyancy_term"))));
            EXPECT_NEAR(
                errors[2], 0.0_rt,
                reference_tol *
                    std::abs(
                        utils::field_max(repo.get_field("sdr_dissipation"))));
            if (srate < 1.0_rt) {
                EXPECT_LT(errors[4], 0.0_rt);
            } else {
                EXPECT_GT(errors[3], 0.0_rt);
            }
            EXPECT_GT(utils::field_max(repo.get_field("shear_prod")), 0.0_rt);
            EXPECT_GT(utils::field_max(repo.get_field("omega_src")), 0.0_rt);
            EXPECT_GT(
                utils::field_max(repo.get_field("sdr_lhs_src_term")), 0.0_rt);
            if (model == "KOmegaSST") {
                EXPECT_GT(
                    std::abs(utils::field_min(repo.get_field("buoyancy_term"))),
                    0.0_rt);
            }
            if (diff_type == DiffusionType::Crank_Nicolson) {
                EXPECT_GT(
                    utils::field_max(repo.get_field("sdr_dissipation")),
                    0.0_rt);
            }

            std::vector<std::unique_ptr<kynema_sgf::ScratchField>> reference;
            for (const auto& entry : fields) {
                reference.emplace_back(repo.create_scratch_field(1, 0));
                kynema_sgf::field_ops::copy(
                    *reference.back(), repo.get_field(entry.first), 0, 0, 1, 0);
            }

            init_linear_z(rho, 2.5_rt, 0.0005_rt);
            init_linear_z(rho_nph, 2.5_rt, 0.0005_rt);
            tmodel.update_turbulent_viscosity(
                kynema_sgf::FieldState::New, diff_type);
            tke_eqn.compute_source_term(kynema_sgf::FieldState::NPH);
            sdr_eqn.compute_source_term(kynema_sgf::FieldState::NPH);

            auto error = repo.create_scratch_field(1, 0);
            for (std::size_t n = 0; n < fields.size(); ++n) {
                const auto& [name, scale] = fields[n];
                SCOPED_TRACE(name);
                kynema_sgf::field_ops::lincomb(
                    *error, 1.0_rt, repo.get_field(name), 0, -scale,
                    *reference[n], 0, 0, 1, 0);
                const amrex::Real magnitude = amrex::max<amrex::Real>(
                    std::abs(utils::field_min(*reference[n])),
                    std::abs(utils::field_max(*reference[n])));
                // The cross-diffusion floor also changes the blending
                // function slightly when density is scaled in single precision.
                const amrex::Real tol =
                    20.0_rt * kynema_sgf::constants::TIGHT_TOL *
                    amrex::max<amrex::Real>(magnitude, 1.0e-6_rt);
                EXPECT_NEAR(utils::field_min(*error), 0.0_rt, tol);
                EXPECT_NEAR(utils::field_max(*error), 0.0_rt, tol);
            }
        }
    }
}

TEST_F(TurbRANSTest, test_kw_sst_density_sources)
{
    check_kw_density_sources("KOmegaSST");
}

TEST_F(TurbRANSTest, test_kw_sst_iddes_density_sources)
{
    check_kw_density_sources("KOmegaSSTIDDES");
}

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
    auto& shear_prod = sim().repo().get_field("shear_prod");

    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e3_rt;
    tke.setVal(tke_val);
    tlscale.setVal(1.0_rt);
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    auto err = max_production_errors(
        shear_prod, Cmu, lambda, kappa, tke_val, rho0, srate, nsq);
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
    err = max_production_errors(
        shear_prod, Cmu, lambda, kappa, tke_val, rho0, srate, nsq);
    EXPECT_NEAR(err[0], 0.0_rt, tol);
    EXPECT_NEAR(err[1], 0.0_rt, tol);
    EXPECT_NEAR(err[2], 0.0_rt, tol);
}

} // namespace kynema_sgf_tests
