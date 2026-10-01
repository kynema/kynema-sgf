#include <numbers>
#include "gtest/gtest.h"
#include "ks_test_utils/MeshTest.H"
#include "src/turbulence/TurbulenceModel.H"
#include "ks_test_utils/test_utils.H"
#include "src/utilities/math_ops.H"
#include "src/utilities/tagging/CartBoxRefinement.H"
#include "src/turbulence/LES/hybrid_length_scale.H"
#include "src/wind_energy/ABL.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParReduce.H"
#include "AMReX_iMultiFab.H"
#include <sstream>

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

void init_field3(kynema_sgf::Field& fld, amrex::Real srate)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();

    amrex::Real offset = 0.0_rt;
    if (fld.field_location() == kynema_sgf::FieldLoc::CELL) {
        offset = 0.5_rt;
    }

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& dx = mesh.Geom(lev).CellSizeArray();
        const auto& problo = mesh.Geom(lev).ProbLoArray();
        const auto& farrs = fld(lev).arrays();

        amrex::ParallelFor(
            fld(lev), fld.num_grow(), fld.num_comp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                const amrex::IntVect iv(i, j, k);
                const amrex::Real xc = problo[n] + ((iv[n] + offset) * dx[n]);
                farrs[nbx](i, j, k, n) = xc / std::sqrt(6.0_rt) * srate;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void init_field_amd(kynema_sgf::Field& fld, amrex::Real scale)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();

    amrex::Real offset = 0.0_rt;
    if (fld.field_location() == kynema_sgf::FieldLoc::CELL) {
        offset = 0.5_rt;
    }

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

                farrs[nbx](i, j, k, 0) = 1 * x / std::sqrt(6.0_rt) * scale;
                farrs[nbx](i, j, k, 1) = -2 * y / std::sqrt(6.0_rt) * scale;
                farrs[nbx](i, j, k, 2) = -1 * z / std::sqrt(6.0_rt) * scale;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void init_field_incomp(kynema_sgf::Field& fld, amrex::Real scale)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();

    amrex::Real offset = 0.0_rt;
    if (fld.field_location() == kynema_sgf::FieldLoc::CELL) {
        offset = 0.5_rt;
    }

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

                farrs[nbx](i, j, k, 0) = 1.0_rt * x * scale;
                farrs[nbx](i, j, k, 1) = -2.0_rt * y * scale;
                farrs[nbx](i, j, k, 2) = 1.0_rt * z * scale;
            });
    }
    amrex::Gpu::streamSynchronize();
}

void init_field1(kynema_sgf::Field& fld, amrex::Real tgrad)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();

    amrex::Real offset = 0.0_rt;
    if (fld.field_location() == kynema_sgf::FieldLoc::CELL) {
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

class TurbLESTest : public MeshTest
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
        }
    }

    const amrex::Real m_dx = 10.0_rt / 10.0_rt;
    const amrex::Real m_dy = 10.0_rt / 20.0_rt;
    const amrex::Real m_dz = 10.0_rt / 30.0_rt;
};

TEST_F(TurbLESTest, test_smag_setup_calc)
{
    // Parser inputs for turbulence model
    const amrex::Real Cs = 0.16_rt;
    const amrex::Real visc = 1.0e-5_rt;
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "Smagorinsky");
    }
    {
        amrex::ParmParse pp("Smagorinsky_coeffs");
        pp.add("Cs", Cs);
    }
    {
        amrex::ParmParse pp("transport");
        pp.add("viscosity", visc);
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

    for (const std::pair<const std::string, const amrex::Real> n : model_dict) {
        // Only a single model parameter, Cs
        EXPECT_EQ(n.first, "Cs");
        EXPECT_EQ(n.second, Cs);
    }

    // Constants for fields
    const amrex::Real srate = 0.5_rt;
    const amrex::Real rho0 = 1.2_rt;

    // Set up velocity field with constant strainrate
    auto& vel = sim().repo().get_field("velocity");
    init_field3(vel, srate);
    // Set up uniform unity density field
    auto& dens = sim().repo().get_field("density");
    dens.setVal(rho0);

    // Update turbulent viscosity directly
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const auto& muturb = sim().repo().get_field("mu_turb");

    // Check values of turbulent viscosity
    auto min_val = utils::field_min(muturb);
    auto max_val = utils::field_max(muturb);
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    const amrex::Real smag_answer =
        rho0 * kynema_sgf::utils::powi(Cs, 2) *
        kynema_sgf::utils::powi(std::cbrt(m_dx * m_dy * m_dz), 2) * srate;
    EXPECT_NEAR(min_val, smag_answer, tol);
    EXPECT_NEAR(max_val, smag_answer, tol);

    // Check values of effective viscosity
    auto& mueff = sim().repo().get_field("velocity_mueff");
    tmodel.update_mueff(mueff);
    min_val = utils::field_min(mueff);
    max_val = utils::field_max(mueff);
    EXPECT_NEAR(min_val, smag_answer + 1.0e-5_rt, tol);
    EXPECT_NEAR(max_val, smag_answer + 1.0e-5_rt, tol);

    // Check that this effective viscosity is what gets to icns diffusion
    auto visc_name = pde_mgr.icns().fields().mueff.name();
    EXPECT_EQ(visc_name, "velocity_mueff");
}

TEST_F(TurbLESTest, test_1eqKsgs_setup_calc)
{
    // Parser inputs for turbulence model
    const amrex::Real Ceps = 0.11_rt;
    const amrex::Real Ce = 0.99_rt;
    const amrex::Real Tref = 263.5_rt;
    const amrex::Real gravz = 10.0_rt;
    const amrex::Real rho0 = 1.2_rt;
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "OneEqKsgsM84");
    }
    {
        amrex::ParmParse pp("OneEqKsgsM84_coeffs");
        pp.add("Ceps", Ceps);
        pp.add("Ce", Ce);
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
        pp.add("surface_temp_rate", -0.25_rt);
        amrex::Vector<amrex::Real> t_hts{0.0_rt, 100.0_rt, 400.0_rt};
        pp.addarr("temperature_heights", t_hts);
        amrex::Vector<amrex::Real> t_vals{265.0_rt, 265.0_rt, 268.0_rt};
        pp.addarr("temperature_values", t_vals);
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
    sim().create_transport_model();
    sim().init_physics();
    sim().create_turbulence_model();
    sim().turbulence_model().post_init_actions();
    auto& tmodel = sim().turbulence_model();

    // Get coefficients
    auto model_dict = tmodel.model_coeffs();

    int ct = 0;
    for (const std::pair<const std::string, const amrex::Real> n : model_dict) {
        // Two model parameters
        if (ct == 0) {
            EXPECT_EQ(n.first, "Ceps");
            EXPECT_EQ(n.second, Ceps);
        } else {
            EXPECT_EQ(n.first, "Ce");
            EXPECT_EQ(n.second, Ce);
        }
        ++ct;
    }

    // Constants for fields
    const amrex::Real srate = 20.0_rt;
    const amrex::Real Tgz = 2.0_rt;
    const amrex::Real tlscale_val = 1.1_rt;
    const amrex::Real tke_val = 0.1_rt;
    // Set up velocity field with constant strainrate
    auto& vel = sim().repo().get_field("velocity");
    init_field3(vel, srate);
    // Set up uniform unity density field
    auto& dens = sim().repo().get_field("density");
    dens.setVal(rho0);
    // Set up temperature field with constant gradient in z
    auto& temp = sim().repo().get_field("temperature");
    init_field1(temp, Tgz);
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
    const auto min_val = utils::field_min(muturb);
    const auto max_val = utils::field_max(muturb);
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    const amrex::Real ksgs_answer =
        rho0 * Ce *
        amrex::min<amrex::Real>(
            std::cbrt(m_dx * m_dy * m_dz),
            0.76_rt * std::sqrt(tke_val / (Tgz * gravz) * Tref)) *
        std::sqrt(tke_val);
    EXPECT_NEAR(min_val, ksgs_answer, tol);
    EXPECT_NEAR(max_val, ksgs_answer, tol);
}

//! Two-level mesh with a refined box that touches the lower wall, as in
//! https://github.com/kynema/kynema-sgf/issues/1806
class TurbLESLevelTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            const amrex::Vector<int> ncell{{m_nx, m_nx, m_nx}};
            pp.addarr("n_cell", ncell);
            pp.add("max_level", 1);
            pp.add("max_grid_size", m_nx);
            pp.add("blocking_factor", 2);
        }
        {
            amrex::ParmParse pp("geometry");
            const amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            const amrex::Vector<amrex::Real> probhi{{m_len, m_len, m_len}};
            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
            pp.addarr("is_periodic", m_periodic);
        }

        // Refine a box in the middle of the domain, from the wall up
        std::stringstream ss;
        ss << "1 // Number of levels" << '\n';
        ss << "1 // Number of boxes at this level" << '\n';
        ss << "16.1 16.1 0.0 47.9 47.9 31.9" << '\n';

        create_mesh_instance<RefineMesh>();
        std::unique_ptr<kynema_sgf::CartBoxRefinement> box_refine(
            new kynema_sgf::CartBoxRefinement(sim()));
        box_refine->read_inputs(mesh(), ss);

        if (mesh<RefineMesh>() != nullptr) {
            mesh<RefineMesh>()->refine_criteria_vec().push_back(
                std::move(box_refine));
        }
    }

    const int m_nx{8};
    const amrex::Real m_len{64.0_rt};
    amrex::Vector<int> m_periodic{{1, 1, 1}};

public:
    //! Filter width cbrt(dx dy dz) on a level
    [[nodiscard]] amrex::Real filter_width(const int lev) const
    {
        return m_len / static_cast<amrex::Real>(m_nx * (1 << lev));
    }
};

// Documents the cause of issue 1806: for the same subgrid kinetic energy the
// OneEqKsgsM84 eddy viscosity halves with every refinement level, because its
// length scale is the local filter width
TEST_F(TurbLESLevelTest, test_1eqKsgs_level_dependence)
{
    const amrex::Real Ceps = 0.93_rt;
    const amrex::Real Ce = 0.1_rt;
    const amrex::Real Tref = 300.0_rt;
    const amrex::Real rho0 = 1.2_rt;
    const amrex::Real tke_val = 0.3_rt;
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "OneEqKsgsM84");
    }
    {
        amrex::ParmParse pp("OneEqKsgsM84_coeffs");
        pp.add("Ceps", Ceps);
        pp.add("Ce", Ce);
    }
    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<std::string> physics{"ABL"};
        pp.addarr("physics", physics);
        pp.add("density", rho0);
        amrex::Vector<amrex::Real> vvec{8.0_rt, 0.0_rt, 0.0_rt};
        pp.addarr("velocity", vvec);
        amrex::Vector<amrex::Real> gvec{0.0_rt, 0.0_rt, -9.81_rt};
        pp.addarr("gravity", gvec);
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("surface_temp_flux", 0.0_rt);
        amrex::Vector<amrex::Real> t_hts{0.0_rt, 100.0_rt};
        pp.addarr("temperature_heights", t_hts);
        amrex::Vector<amrex::Real> t_vals{Tref, Tref};
        pp.addarr("temperature_values", t_vals);
    }
    {
        amrex::ParmParse pp("transport");
        pp.add("reference_temperature", Tref);
    }

    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().create_transport_model();
    sim().init_physics();
    sim().create_turbulence_model();
    sim().turbulence_model().post_init_actions();
    auto& tmodel = sim().turbulence_model();

    ASSERT_EQ(sim().repo().num_active_levels(), 2);

    // Neutral, uniform state: the same subgrid kinetic energy on both levels
    sim().repo().get_field("velocity").setVal(0.0_rt);
    sim().repo().get_field("density").setVal(rho0);
    sim().repo().get_field("temperature").setVal(Tref);
    sim().repo().get_field("tke").setVal(tke_val);

    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const auto& muturb = sim().repo().get_field("mu_turb");

    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    for (int lev = 0; lev < 2; ++lev) {
        const amrex::Real answer =
            rho0 * Ce * filter_width(lev) * std::sqrt(tke_val);
        EXPECT_NEAR(muturb(lev).min(0), answer, tol);
        EXPECT_NEAR(muturb(lev).max(0), answer, tol);
    }
    // Same physical location, same subgrid kinetic energy, half the eddy
    // viscosity on the finer level
    EXPECT_NEAR(muturb(1).max(0) / muturb(0).max(0), 0.5_rt, tol);
}

// Shows that the near-wall hybrid RANS-LES length scale shared with the
// Kosovic model removes most of that level dependence close to the wall when
// it is evaluated with the OneEqKsgsM84 coefficients
TEST_F(TurbLESLevelTest, test_hybrid_length_level_independence)
{
    namespace hl = kynema_sgf::turbulence::hybrid_length;
    const amrex::Real Ceps = 0.93_rt;
    const amrex::Real Ce = 0.1_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real switch_height = 24.0_rt;
    const amrex::Real exponent = 2.0_rt;
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    populate_parameters();
    initialize_mesh();
    ASSERT_EQ(sim().repo().num_active_levels(), 2);
    const amrex::Real ds0 = filter_width(0);
    const amrex::Real ds1 = filter_width(1);

    // Limits of the blending
    EXPECT_NEAR(hl::rans_weight(0.0_rt, switch_height), 1.0_rt, tol);
    EXPECT_NEAR(
        hl::blended_length_sqr(ds0 * ds0, 2.0_rt, 0.0_rt, exponent), ds0 * ds0,
        tol);
    EXPECT_NEAR(
        hl::blended_length_sqr(ds0 * ds0, 2.0_rt, 1.0_rt, exponent), 2.0_rt,
        tol);

    // The RANS length scale recovers the log-law eddy viscosity when
    // production balances dissipation: k = (Ce / Ceps) l^2 S^2
    {
        const amrex::Real utau = 0.4_rt;
        const amrex::Real z = 10.0_rt;
        const amrex::Real l_rans =
            hl::one_eq_rans_length(kappa, z, 1.0_rt, Ce, Ceps);
        const amrex::Real shear = utau / (kappa * z);
        const amrex::Real k_eq = (Ce / Ceps) * l_rans * l_rans * shear * shear;
        EXPECT_NEAR(Ce * l_rans * std::sqrt(k_eq), kappa * utau * z, tol);
    }

    // Ratio of the fine to the coarse length scale (and so of the eddy
    // viscosity for the same subgrid kinetic energy) at the same height
    const auto level_ratio = [&](const amrex::Real z, const amrex::Real ds_c,
                                 const amrex::Real ds_f) {
        const amrex::Real l_rans =
            hl::one_eq_rans_length(kappa, z, 1.0_rt, Ce, Ceps);
        const amrex::Real w = hl::rans_weight(z, switch_height);
        const amrex::Real l_c = std::sqrt(
            hl::blended_length_sqr(ds_c * ds_c, l_rans * l_rans, w, exponent));
        const amrex::Real l_f = std::sqrt(
            hl::blended_length_sqr(ds_f * ds_f, l_rans * l_rans, w, exponent));
        return l_f / l_c;
    };

    // Close to the wall the two levels agree (OneEqKsgsM84 gives 0.5)
    EXPECT_GT(level_ratio(0.5_rt * ds0, ds0, ds1), 0.95_rt);
    EXPECT_GT(level_ratio(2.0_rt * ds0, ds0, ds1), 0.9_rt);
    // Far above the wall the LES filter width is recovered on each level
    EXPECT_NEAR(level_ratio(20.0_rt * switch_height, ds0, ds1), ds1 / ds0, tol);

    // Filter widths of a 64 x 64 x 16 m grid and of its second refinement
    // level: OneEqKsgsM84 gives 0.25 at every height
    const amrex::Real ds_coarse = std::cbrt(64.0_rt * 64.0_rt * 16.0_rt);
    const amrex::Real ds_fine = 0.25_rt * ds_coarse;
    EXPECT_GT(level_ratio(8.0_rt, ds_coarse, ds_fine), 0.75_rt);
    EXPECT_GT(level_ratio(16.0_rt, ds_coarse, ds_fine), 0.7_rt);
}

namespace {

//! Fill the velocity with a logarithmic wind along x that depends on the
//! height only. The ghost cells below the wall repeat the first cell so that
//! every plane average reads finite values.
void init_log_wind(
    kynema_sgf::Field& vel,
    const amrex::Real utau,
    const amrex::Real kappa,
    const amrex::Real z0)
{
    const auto& mesh = vel.repo().mesh();
    const int nlevels = vel.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const amrex::Real dz = mesh.Geom(lev).CellSizeArray()[2];
        const amrex::Real zlo = mesh.Geom(lev).ProbLoArray()[2];
        const auto& varrs = vel(lev).arrays();
        amrex::ParallelFor(
            vel(lev), vel.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const int kk = amrex::max(k, 0);
                const amrex::Real z = zlo + ((kk + 0.5_rt) * dz);
                varrs[nbx](i, j, k, 0) = utau / kappa * std::log(z / z0);
                varrs[nbx](i, j, k, 1) = 0.0_rt;
                varrs[nbx](i, j, k, 2) = 0.0_rt;
            });
    }
    amrex::Gpu::streamSynchronize();
}

//! Cells of a level not covered by a finer level (1) or covered (0)
amrex::iMultiFab level_mask(const kynema_sgf::FieldRepo& repo, const int lev)
{
    const auto& mesh = repo.mesh();
    amrex::iMultiFab mask;
    if (lev < repo.num_active_levels() - 1) {
        mask = amrex::makeFineMask(
            mesh.boxArray(lev), mesh.DistributionMap(lev),
            mesh.boxArray(lev + 1), mesh.refRatio(lev), 1, 0);
    } else {
        mask.define(mesh.boxArray(lev), mesh.DistributionMap(lev), 1, 0);
        mask.setVal(1);
    }
    return mask;
}

//! Mean over the uncovered wall cells of a level of the x-stress the wall
//! model applies: the ghost cell below the wall holds the wall-normal
//! gradient, so the stress is ghost * mueff / rho of the wall cell
amrex::Real wall_stress_mean(
    const kynema_sgf::Field& vel,
    const kynema_sgf::Field& mueff,
    const kynema_sgf::Field& rho,
    const int lev,
    amrex::Real& ncells)
{
    const auto& repo = vel.repo();
    const int kwall = repo.mesh().Geom(lev).Domain().smallEnd(2);
    const amrex::iMultiFab mask = level_mask(repo, lev);
    const auto& v_arrs = vel(lev).const_arrays();
    const auto& mu_arrs = mueff(lev).const_arrays();
    const auto& rho_arrs = rho(lev).const_arrays();
    const auto& msk_arrs = mask.const_arrays();

    using SumTuple = amrex::GpuTuple<amrex::Real, amrex::Real>;
    const SumTuple sums = amrex::ParReduce(
        amrex::TypeList<amrex::ReduceOpSum, amrex::ReduceOpSum>{},
        amrex::TypeList<amrex::Real, amrex::Real>{}, vel(lev),
        amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) -> SumTuple {
            if (k != kwall) {
                return {0.0_rt, 0.0_rt};
            }
            const auto msk = static_cast<amrex::Real>(msk_arrs[nbx](i, j, k));
            const amrex::Real tau = v_arrs[nbx](i, j, k - 1, 0) *
                                    mu_arrs[nbx](i, j, k) /
                                    rho_arrs[nbx](i, j, k);
            return {msk * tau, msk};
        });
    amrex::GpuArray<amrex::Real, 2> vals{
        amrex::get<0>(sums), amrex::get<1>(sums)};
    amrex::ParallelDescriptor::ReduceRealSum(vals.data(), 2);
    ncells = vals[1];
    return vals[0] / amrex::max(vals[1], 1.0_rt);
}

//! Mean over the uncovered cells (i, j, kcell) of a level of the modeled
//! x-stress across the face above the cell, mu_t du/dz / rho, with the
//! viscosity averaged from the two cells that share the face
amrex::Real sgs_stress_mean(
    const kynema_sgf::Field& vel,
    const kynema_sgf::Field& mu_turb,
    const kynema_sgf::Field& rho,
    const int lev,
    const int kcell)
{
    const auto& repo = vel.repo();
    const amrex::Real dz = repo.mesh().Geom(lev).CellSizeArray()[2];
    const amrex::iMultiFab mask = level_mask(repo, lev);
    const auto& v_arrs = vel(lev).const_arrays();
    const auto& mu_arrs = mu_turb(lev).const_arrays();
    const auto& rho_arrs = rho(lev).const_arrays();
    const auto& msk_arrs = mask.const_arrays();

    using SumTuple = amrex::GpuTuple<amrex::Real, amrex::Real>;
    const SumTuple sums = amrex::ParReduce(
        amrex::TypeList<amrex::ReduceOpSum, amrex::ReduceOpSum>{},
        amrex::TypeList<amrex::Real, amrex::Real>{}, vel(lev),
        amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) -> SumTuple {
            if (k != kcell) {
                return {0.0_rt, 0.0_rt};
            }
            const auto msk = static_cast<amrex::Real>(msk_arrs[nbx](i, j, k));
            const amrex::Real mu_face =
                0.5_rt * (mu_arrs[nbx](i, j, k) + mu_arrs[nbx](i, j, k + 1));
            const amrex::Real dudz =
                (v_arrs[nbx](i, j, k + 1, 0) - v_arrs[nbx](i, j, k, 0)) / dz;
            return {msk * mu_face * dudz / rho_arrs[nbx](i, j, k), msk};
        });
    amrex::GpuArray<amrex::Real, 2> vals{
        amrex::get<0>(sums), amrex::get<1>(sums)};
    amrex::ParallelDescriptor::ReduceRealSum(vals.data(), 2);
    return vals[0] / amrex::max(vals[1], 1.0_rt);
}

} // namespace

/** Two-level mesh of TurbLESLevelTest with a wall-modeled lower boundary */
class TurbLESWallLevelTest : public TurbLESLevelTest
{
protected:
    void populate_parameters() override
    {
        m_periodic = {1, 1, 0};
        TurbLESLevelTest::populate_parameters();
        {
            amrex::ParmParse pp("zlo");
            pp.add("type", (std::string) "wall_model");
            pp.add("temperature_type", (std::string) "wall_model");
        }
        {
            amrex::ParmParse pp("zhi");
            pp.add("type", (std::string) "slip_wall");
            pp.add("temperature_type", (std::string) "fixed_gradient");
        }
    }

public:
    //! Stress budget of the first coarse cell inside and outside the box,
    //! with the filter width (surface_rans = false) or the near-wall hybrid
    //! length scale (surface_rans = true) of OneEqKsgsM84
    void run_stress_budget(bool surface_rans);
};

/** Stress budget of the first coarse cell inside and outside a refined box
 *
 *  A logarithmic wind with friction velocity utau and the subgrid kinetic
 *  energy of the equilibrium log layer are set on both levels. The wall model
 *  (local stress, whose magnitude is the plane-mean utau^2) then applies the
 *  same stress under the coarse and the fine wall cells. Across the plane at
 *  the top of the coarse wall cell, z = dz0, the modeled stress of
 *  OneEqKsgsM84 in the refined box is less than half of the coarse value for
 *  the same velocity and the same subgrid kinetic energy, because its length
 *  scale is the filter width of each level. The wall-distance length scale of
 *  the hybrid helper gives kappa utau z on both levels, so the modeled stress
 *  at z = dz0 is the wall stress on both, up to the discrete log law.
 */
void TurbLESWallLevelTest::run_stress_budget(const bool surface_rans)
{
    namespace hl = kynema_sgf::turbulence::hybrid_length;
    const amrex::Real Ceps = 0.93_rt;
    const amrex::Real Ce = 0.1_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real z0 = 0.1_rt;
    const amrex::Real utau = 0.4_rt;
    const amrex::Real utau2 = utau * utau;
    const amrex::Real Tref = 300.0_rt;
    const amrex::Real rho0 = 1.2_rt;
    const amrex::Real mu = 0.01_rt;
    // Subgrid kinetic energy of the log layer when production balances
    // dissipation: nu_t S = utau^2, S = utau / (kappa z) and
    // nu_t S^2 = Ceps k^(3/2) / l with nu_t = Ce l sqrt(k) give
    // k = utau^2 / sqrt(Ce Ceps)
    const amrex::Real tke_eq = utau2 / std::sqrt(Ce * Ceps);
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    // Log wind at the reference height of the wall function, 0.5 dz0
    const amrex::Real u_ref =
        utau / kappa * std::log(0.5_rt * filter_width(0) / z0);
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "OneEqKsgsM84");
    }
    {
        amrex::ParmParse pp("OneEqKsgsM84_coeffs");
        pp.add("Ceps", Ceps);
        pp.add("Ce", Ce);
    }
    if (surface_rans) {
        // Near-wall hybrid length scale, switch height of three coarse cells
        amrex::ParmParse pp("OneEqKsgsM84");
        pp.add("surfaceRANS", 1);
        pp.add("switchLoc", 3.0_rt * filter_width(0));
        pp.add("surfaceRANSExp", 2.0_rt);
    }
    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<std::string> physics{"ABL"};
        pp.addarr("physics", physics);
        pp.add("density", rho0);
        amrex::Vector<amrex::Real> vvec{8.0_rt, 0.0_rt, 0.0_rt};
        pp.addarr("velocity", vvec);
        amrex::Vector<amrex::Real> gvec{0.0_rt, 0.0_rt, -9.81_rt};
        pp.addarr("gravity", gvec);
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("kappa", kappa);
        pp.add("surface_roughness_z0", z0);
        pp.add("surface_temp_flux", 0.0_rt);
        pp.add("wall_shear_stress_type", (std::string) "local");
        // Hand the wall function its reference velocity instead of the
        // plane average, so that the test does not depend on the averaging
        // stencil next to the wall of a mesh refined at the wall
        pp.add("inflow_outflow_mode", 1);
        amrex::Vector<amrex::Real> wf_vel{u_ref, 0.0_rt};
        pp.addarr("wf_velocity", wf_vel);
        pp.add("wf_vmag", u_ref);
        pp.add("wf_theta", Tref);
        amrex::Vector<amrex::Real> t_hts{0.0_rt, 100.0_rt};
        pp.addarr("temperature_heights", t_hts);
        amrex::Vector<amrex::Real> t_vals{Tref, Tref};
        pp.addarr("temperature_values", t_vals);
    }
    {
        amrex::ParmParse pp("transport");
        pp.add("reference_temperature", Tref);
        pp.add("viscosity", mu);
    }

    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().create_transport_model();
    sim().init_physics();
    sim().create_turbulence_model();
    sim().turbulence_model().post_init_actions();
    auto& tmodel = sim().turbulence_model();

    ASSERT_EQ(sim().repo().num_active_levels(), 2);
    const amrex::Real dz0 = filter_width(0);
    const amrex::Real dz1 = filter_width(1);

    auto& repo = sim().repo();
    auto& velocity = repo.get_field("velocity");
    auto& density = repo.get_field("density");
    auto& vel_mueff = repo.get_field("velocity_mueff");
    init_log_wind(velocity, utau, kappa, z0);
    density.setVal(rho0);
    repo.get_field("temperature").setVal(Tref);
    repo.get_field("tke").setVal(tke_eq);
    vel_mueff.setVal(mu);

    // Wall function: the log wind at zref = 0.5 dz0 gives back the friction
    // velocity of the profile
    for (auto& phys : sim().physics()) {
        phys->post_init_actions();
    }
    const auto& abl = sim().physics_manager().get<kynema_sgf::ABL>();
    const auto& mo = abl.abl_wall_function().mo();
    EXPECT_NEAR(mo.zref, 0.5_rt * dz0, tol);
    EXPECT_NEAR(mo.utau, utau, tol);

    // Fill the wall-model ghost cells on both levels
    pde_mgr.advance_states();
    velocity.apply_bc_funcs(kynema_sgf::FieldState::Old);

    // The wall model applies the same stress utau^2 under the coarse and the
    // fine wall cells
    for (int lev = 0; lev < 2; ++lev) {
        amrex::Real ncells = 0.0_rt;
        const amrex::Real tau_w =
            wall_stress_mean(velocity, vel_mueff, density, lev, ncells);
        ASSERT_GT(ncells, 0.0_rt) << "level " << lev;
        EXPECT_NEAR(tau_w, utau2, tol * utau2) << "level " << lev;
    }

    // OneEqKsgsM84: nu_t = Ce ds sqrt(k) with the filter width of each level
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const auto& mu_turb = repo.get_field("mu_turb");
    if (!surface_rans) {
        for (int lev = 0; lev < 2; ++lev) {
            const amrex::Real answer =
                rho0 * Ce * filter_width(lev) * std::sqrt(tke_eq);
            EXPECT_NEAR(mu_turb(lev).min(0), answer, tol);
            EXPECT_NEAR(mu_turb(lev).max(0), answer, tol);
        }
    }

    // Discrete shear of the log wind across the plane z = dz0: between the
    // level-0 cells centered at 0.5 dz0 and 1.5 dz0, and between the level-1
    // cells centered at 1.5 dz1 and 2.5 dz1
    const amrex::Real shear0 = utau / kappa * std::log(3.0_rt) / dz0;
    const amrex::Real shear1 = utau / kappa * std::log(5.0_rt / 3.0_rt) / dz1;

    // Modeled stress across z = dz0: the level-0 wall cell (kcell = 0) and
    // the second level-1 cell (kcell = 1) share that plane
    const amrex::Real sgs0 = sgs_stress_mean(velocity, mu_turb, density, 0, 0);
    const amrex::Real sgs1 = sgs_stress_mean(velocity, mu_turb, density, 1, 1);
    if (surface_rans) {
        // With the hybrid length scale (switch height 3 dz0) the modeled
        // stress at z = dz0 is about 0.7 u_*^2 on both levels: the box is no
        // longer under-mixed relative to its surroundings
        EXPECT_GT(sgs0, 0.6_rt * utau2);
        EXPECT_GT(sgs1, 0.6_rt * utau2);
        EXPECT_GT(sgs1 / sgs0, 0.9_rt);
        EXPECT_LT(sgs1 / sgs0, 1.1_rt);
    } else {
        EXPECT_NEAR(sgs0, Ce * dz0 * std::sqrt(tke_eq) * shear0, tol);
        EXPECT_NEAR(sgs1, Ce * dz1 * std::sqrt(tke_eq) * shear1, tol);
        // Same velocity, same subgrid kinetic energy, same wall stress
        // below: the refined box transmits ln(5/3) / ln(3) = 0.47 of the
        // coarse stress
        EXPECT_NEAR(
            sgs1 / sgs0, std::log(5.0_rt / 3.0_rt) / std::log(3.0_rt), tol);
        EXPECT_LT(sgs1 / sgs0, 0.5_rt);
    }

    // With the wall-distance length scale at the height of that plane the
    // eddy viscosity is kappa utau z on both levels, and the modeled stress
    // is the wall stress on both up to the discrete log law (within 10 %,
    // fine to coarse ratio 2 ln(5/3) / ln(3) = 0.93)
    const amrex::Real nu_rans =
        Ce * hl::one_eq_rans_length(kappa, dz0, 1.0_rt, Ce, Ceps) *
        std::sqrt(tke_eq);
    EXPECT_NEAR(nu_rans, kappa * utau * dz0, tol);
    EXPECT_NEAR(nu_rans * shear0, utau2, 0.1_rt * utau2);
    EXPECT_NEAR(nu_rans * shear1, utau2, 0.1_rt * utau2);
    EXPECT_GT(shear1 / shear0, 0.9_rt);
}

TEST_F(TurbLESWallLevelTest, test_1eqKsgs_wall_stress_budget)
{
    run_stress_budget(false);
}

TEST_F(TurbLESWallLevelTest, test_1eqKsgs_wall_stress_budget_surface_rans)
{
    run_stress_budget(true);
}

namespace {

//! Value of a cell-centered field in one cell
amrex::Real cell_value(
    const kynema_sgf::Field& fld,
    const int lev,
    const int i,
    const int j,
    const int k)
{
    const amrex::IntVect iv(i, j, k);
    return fld(lev).max(amrex::Box(iv, iv), 0);
}

} // namespace

/** OneEqKsgsM84.surfaceRANS: the length scale, the eddy viscosity, the
 *  diffusivity and the TKE dissipation follow the hybrid blend on both
 *  levels, and the box is not under-mixed relative to the coarse level
 */
TEST_F(TurbLESLevelTest, test_1eqKsgs_surface_rans_level_independence)
{
    namespace hl = kynema_sgf::turbulence::hybrid_length;
    const amrex::Real Ceps = 0.93_rt;
    const amrex::Real Ce = 0.1_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real switch_loc = 24.0_rt;
    const amrex::Real exponent = 2.0_rt;
    const amrex::Real Tref = 300.0_rt;
    const amrex::Real rho0 = 1.2_rt;
    const amrex::Real tke_val = 0.3_rt;
    const amrex::Real mu = 1.0e-5_rt;
    const amrex::Real prandtl = 0.7_rt;
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "OneEqKsgsM84");
    }
    {
        amrex::ParmParse pp("OneEqKsgsM84_coeffs");
        pp.add("Ceps", Ceps);
        pp.add("Ce", Ce);
    }
    {
        amrex::ParmParse pp("OneEqKsgsM84");
        pp.add("surfaceRANS", 1);
        pp.add("switchLoc", switch_loc);
        pp.add("surfaceRANSExp", exponent);
    }
    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<std::string> physics{"ABL"};
        pp.addarr("physics", physics);
        pp.add("density", rho0);
        amrex::Vector<amrex::Real> vvec{8.0_rt, 0.0_rt, 0.0_rt};
        pp.addarr("velocity", vvec);
        amrex::Vector<amrex::Real> gvec{0.0_rt, 0.0_rt, -9.81_rt};
        pp.addarr("gravity", gvec);
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("kappa", kappa);
        pp.add("surface_temp_flux", 0.0_rt);
        amrex::Vector<amrex::Real> t_hts{0.0_rt, 100.0_rt};
        pp.addarr("temperature_heights", t_hts);
        amrex::Vector<amrex::Real> t_vals{Tref, Tref};
        pp.addarr("temperature_values", t_vals);
    }
    {
        amrex::ParmParse pp("transport");
        pp.add("reference_temperature", Tref);
        pp.add("viscosity", mu);
        pp.add("laminar_prandtl", prandtl);
    }

    populate_parameters();
    initialize_mesh();
    auto& pde_mgr = sim().pde_manager();
    pde_mgr.register_icns();
    sim().create_transport_model();
    sim().init_physics();
    sim().create_turbulence_model();
    sim().turbulence_model().post_init_actions();
    auto& tmodel = sim().turbulence_model();
    auto& repo = sim().repo();
    ASSERT_EQ(repo.num_active_levels(), 2);

    // Neutral, uniform state: the same subgrid kinetic energy on both levels
    repo.get_field("velocity").setVal(0.0_rt);
    repo.get_field("density").setVal(rho0);
    repo.get_field("temperature").setVal(Tref);
    repo.get_field("tke").setVal(tke_val);

    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    tmodel.update_alphaeff(repo.get_field("temperature_mueff"));
    for (auto& eqn : pde_mgr.scalar_eqns()) {
        if (eqn->fields().field.name() == "tke") {
            eqn->initialize();
            eqn->compute_source_term(kynema_sgf::FieldState::New);
        }
    }

    const auto& lscale = repo.get_field("hybrid_length_scale");
    const auto& tlscale = repo.get_field("turb_lscale");
    const auto& muturb = repo.get_field("mu_turb");
    const auto& alphaeff = repo.get_field("temperature_mueff");
    const auto& dissip = repo.get_field("dissipation");

    // Neutral (no Obukhov length given, wall function neutral): phi_m = 1
    const auto blend = [&](const amrex::Real ds, const amrex::Real z) {
        const amrex::Real l_rans =
            hl::one_eq_rans_length(kappa, z, 1.0_rt, Ce, Ceps);
        const amrex::Real w = hl::rans_weight(z, switch_loc);
        return std::sqrt(
            hl::blended_length_sqr(ds * ds, l_rans * l_rans, w, exponent));
    };
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    // Every cell of both levels follows the blend at its own height; the
    // eddy viscosity, the diffusivity (1 + 2 l/l = 3) and the dissipation
    // (Ce_local = Ceps for l = ds) use that length
    for (int lev = 0; lev < 2; ++lev) {
        const amrex::Real ds = filter_width(lev);
        const int i = (lev == 0) ? 0 : 7;
        for (int k = 0; k < 8; ++k) {
            const amrex::Real z = (k + 0.5_rt) * ds;
            const amrex::Real l = blend(ds, z);
            EXPECT_NEAR(cell_value(lscale, lev, i, i, k), l, tol)
                << "level " << lev << " k " << k;
            EXPECT_NEAR(cell_value(tlscale, lev, i, i, k), l, tol);
            const amrex::Real mut = rho0 * Ce * l * std::sqrt(tke_val);
            EXPECT_NEAR(cell_value(muturb, lev, i, i, k), mut, tol);
            EXPECT_NEAR(
                cell_value(alphaeff, lev, i, i, k), mu / prandtl + 3.0_rt * mut,
                tol);
            EXPECT_NEAR(
                cell_value(dissip, lev, i, i, k),
                Ceps * tke_val * std::sqrt(tke_val) / l, tol);
        }
    }

    // Level independence near the wall: the first cells of both levels sit
    // close to the log-law length, and the fine level at 6 m has a larger
    // length than the coarse level at 4 m (plain OneEqKsgsM84: half)
    const amrex::Real l0 = cell_value(lscale, 0, 0, 0, 0); // z = 4 m
    const amrex::Real l1 = cell_value(lscale, 1, 7, 7, 1); // z = 6 m
    EXPECT_GT(
        l0 / hl::one_eq_rans_length(kappa, 4.0_rt, 1.0_rt, Ce, Ceps), 0.85_rt);
    EXPECT_GT(
        cell_value(lscale, 1, 7, 7, 0) /
            hl::one_eq_rans_length(kappa, 2.0_rt, 1.0_rt, Ce, Ceps),
        0.85_rt);
    EXPECT_GT(l1 / l0, 1.0_rt);
    // Aloft the blend relaxes toward the filter width of each level
    EXPECT_LT(cell_value(lscale, 0, 0, 0, 7) / filter_width(0), 1.7_rt);
    EXPECT_GT(cell_value(lscale, 0, 0, 0, 7) / filter_width(0), 1.0_rt);
}

/** OneEqKsgsM84.surfaceRANS with ABL.monin_obukhov_length: the stability
 *  function of the wall-distance length uses the given Obukhov length
 */
TEST_F(TurbLESLevelTest, test_1eqKsgs_surface_rans_stable)
{
    namespace hl = kynema_sgf::turbulence::hybrid_length;
    const amrex::Real Ceps = 0.93_rt;
    const amrex::Real Ce = 0.1_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real switch_loc = 24.0_rt;
    const amrex::Real obukhov_len = 100.0_rt;
    const amrex::Real Tref = 300.0_rt;
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "OneEqKsgsM84");
    }
    {
        amrex::ParmParse pp("OneEqKsgsM84_coeffs");
        pp.add("Ceps", Ceps);
        pp.add("Ce", Ce);
    }
    {
        amrex::ParmParse pp("OneEqKsgsM84");
        pp.add("surfaceRANS", 1);
        pp.add("switchLoc", switch_loc);
    }
    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<std::string> physics{"ABL"};
        pp.addarr("physics", physics);
        pp.add("density", 1.0_rt);
        amrex::Vector<amrex::Real> vvec{8.0_rt, 0.0_rt, 0.0_rt};
        pp.addarr("velocity", vvec);
        amrex::Vector<amrex::Real> gvec{0.0_rt, 0.0_rt, -9.81_rt};
        pp.addarr("gravity", gvec);
    }
    {
        amrex::ParmParse pp("ABL");
        pp.add("kappa", kappa);
        pp.add("monin_obukhov_length", obukhov_len);
        pp.add("surface_temp_flux", 0.0_rt);
        amrex::Vector<amrex::Real> t_hts{0.0_rt, 100.0_rt};
        pp.addarr("temperature_heights", t_hts);
        amrex::Vector<amrex::Real> t_vals{Tref, Tref};
        pp.addarr("temperature_values", t_vals);
    }
    {
        amrex::ParmParse pp("transport");
        pp.add("reference_temperature", Tref);
    }

    populate_parameters();
    initialize_mesh();
    sim().pde_manager().register_icns();
    sim().create_transport_model();
    sim().init_physics();
    sim().create_turbulence_model();
    sim().turbulence_model().post_init_actions();
    auto& repo = sim().repo();
    ASSERT_EQ(repo.num_active_levels(), 2);
    repo.get_field("velocity").setVal(0.0_rt);
    repo.get_field("density").setVal(1.0_rt);
    repo.get_field("temperature").setVal(Tref);
    repo.get_field("tke").setVal(0.3_rt);
    sim().turbulence_model().update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);

    const auto& lscale = repo.get_field("hybrid_length_scale");
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    for (int lev = 0; lev < 2; ++lev) {
        const amrex::Real ds = filter_width(lev);
        const int i = (lev == 0) ? 0 : 7;
        for (int k = 0; k < 4; ++k) {
            const amrex::Real z = (k + 0.5_rt) * ds;
            // Stable: phi_m = 1 + 5 z / L shortens the wall-distance length
            const amrex::Real phi_m = 1.0_rt + (5.0_rt * z / obukhov_len);
            const amrex::Real l_rans =
                hl::one_eq_rans_length(kappa, z, phi_m, Ce, Ceps);
            const amrex::Real w = hl::rans_weight(z, switch_loc);
            const amrex::Real l = std::sqrt(
                hl::blended_length_sqr(ds * ds, l_rans * l_rans, w, 2.0_rt));
            EXPECT_NEAR(cell_value(lscale, lev, i, i, k), l, tol)
                << "level " << lev << " k " << k;
        }
    }
}

TEST_F(TurbLESTest, test_AMD_setup_calc)
{
    // Parser inputs for turbulence model
    const amrex::Real C = 0.3_rt;
    const amrex::Real Tref = 200.0_rt;
    const amrex::Real gravz = 10.0_rt;
    const amrex::Real rho0 = 1.0_rt;
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "AMD");
    }
    {
        amrex::ParmParse pp("AMD_coeffs");
        pp.add("C_poincare", C);
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
        amrex::Vector<amrex::Real> t_hts{0.0_rt, 100.0_rt, 400.0_rt};
        pp.addarr("temperature_heights", t_hts);
        amrex::Vector<amrex::Real> t_vals{200.0_rt, 200.0_rt, 200.0_rt};
        pp.addarr("temperature_values", t_vals);
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

    for (const std::pair<const std::string, const amrex::Real> n : model_dict) {
        // Only a single model parameter, Cs
        EXPECT_EQ(n.first, "C_poincare");
        EXPECT_EQ(n.second, C);
    }

    // Constants for fields
    const amrex::Real scale = 1.50_rt;
    const amrex::Real Tgz = 20.0_rt;
    // Set up velocity field with constant strainrate
    auto& vel = sim().repo().get_field("velocity");
    init_field_amd(vel, scale);
    // Set up uniform unity density field
    auto& dens = sim().repo().get_field("density");
    dens.setVal(rho0);
    // Set up temperature field with constant gradient in z
    auto& temp = sim().repo().get_field("temperature");
    init_field1(temp, Tgz);

    // Update turbulent viscosity directly
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const auto& muturb = sim().repo().get_field("mu_turb");

    // Check values of turbulent viscosity
    const auto min_val = utils::field_min(muturb);
    const auto max_val = utils::field_max(muturb);
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    const amrex::Real amd_answer =
        C *
        (-1.0_rt * kynema_sgf::utils::powi(scale / std::sqrt(6.0_rt), 3) *
         (m_dx * m_dx - 8.0_rt * m_dy * m_dy - m_dz * m_dz)) /
        (1.0_rt * scale * scale);
    EXPECT_NEAR(min_val, amd_answer, tol);
    EXPECT_NEAR(max_val, amd_answer, tol);

    // Check values of alphaeff
    auto& alphaeff = sim().repo().declare_cc_field("alphaeff");
    tmodel.update_alphaeff(alphaeff);
    const auto ae_min_val = utils::field_min(alphaeff);
    const auto ae_max_val = utils::field_max(alphaeff);
    const amrex::Real amd_ae_answer =
        C * m_dz * m_dz * scale * 1.0_rt / std::sqrt(6.0_rt);
    EXPECT_NEAR(ae_min_val, amd_ae_answer, tol);
    EXPECT_NEAR(ae_max_val, amd_ae_answer, tol);
}

TEST_F(TurbLESTest, test_AMDNoTherm_setup_calc)
{
    // Parser inputs for turbulence model
    const amrex::Real C = 0.3_rt;
    const amrex::Real rho0 = 1.0_rt;
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "AMDNoTherm");
    }
    {
        amrex::ParmParse pp("AMDNoTherm_coeffs");
        pp.add("C_poincare", C);
    }
    {
        amrex::ParmParse pp("incflo");
        pp.add("density", rho0);
        amrex::Vector<amrex::Real> vvec{8.0_rt, 0.0_rt, 0.0_rt};
        pp.addarr("velocity", vvec);
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

    for (const std::pair<const std::string, const amrex::Real> n : model_dict) {
        // Only a single model parameter, Cs
        EXPECT_EQ(n.first, "C_poincare");
        EXPECT_EQ(n.second, C);
    }

    // Constants for fields
    const amrex::Real scale = 2.0_rt;
    // Set up velocity field with constant strainrate
    auto& vel = sim().repo().get_field("velocity");
    init_field_incomp(vel, scale);
    // Set up uniform unity density field
    auto& dens = sim().repo().get_field("density");
    dens.setVal(rho0);

    // Update turbulent viscosity directly
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const auto& muturb = sim().repo().get_field("mu_turb");

    // Check values of turbulent viscosity
    const auto min_val = utils::field_min(muturb);
    const auto max_val = utils::field_max(muturb);
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

    const amrex::Real amd_answer =
        -C * kynema_sgf::utils::powi(scale, 3) *
        (m_dx * m_dx - 8.0_rt * m_dy * m_dy + m_dz * m_dz) /
        (6 * scale * scale);
    EXPECT_NEAR(min_val, amd_answer, tol);
    EXPECT_NEAR(max_val, amd_answer, tol);
}
TEST_F(TurbLESTest, test_kosovic_setup_calc)
{
    // Parser inputs for turbulence model
    const amrex::Real Cb = 0.36_rt;
    const amrex::Real visc = 1.0e-5_rt;
    const amrex::Real kosovic_Cs = std::sqrt(
        8.0_rt * (1.0_rt + Cb) /
        (27.0_rt * std::numbers::pi_v<amrex::Real> *
         std::numbers::pi_v<amrex::Real>));
    {
        amrex::ParmParse pp("turbulence");
        pp.add("model", (std::string) "Kosovic");
    }
    {
        amrex::ParmParse pp("Kosovic_coeffs");
        pp.add("Cb", Cb);
    }
    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<std::string> physics{"ABL"};
        pp.addarr("physics", physics);
        pp.add("density", 1.225);
        amrex::Vector<amrex::Real> vvec{8.0, 0.0, 0.0};
        pp.addarr("velocity", vvec);
        amrex::Vector<amrex::Real> gvec{0.0, 0.0, -9.81};
        pp.addarr("gravity", gvec);
    }
    {
        amrex::ParmParse pp("transport");
        pp.add("viscosity", visc);
    }
    {
        amrex::ParmParse pp("ABL");
        amrex::Vector<amrex::Real> t_hts{0.0, 100.0, 4000.0};
        pp.addarr("temperature_heights", t_hts);
        amrex::Vector<amrex::Real> t_vals{300.0, 300.0, 300.0};
        pp.addarr("temperature_values", t_vals);
    }
    // Transport
    {
        amrex::ParmParse pp("transport");
        pp.add("reference_temperature", 300.0);
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

    for (const std::pair<const std::string, const amrex::Real> n : model_dict) {
        // Only a single model parameter, Cb
        EXPECT_EQ(n.first, "Cb");
        EXPECT_EQ(n.second, Cb);
    }

    // Constants for fields
    const amrex::Real srate = 0.5_rt;
    const amrex::Real rho0 = 1.2_rt;

    // Set up velocity field with constant strainrate
    auto& vel = sim().repo().get_field("velocity");
    init_field3(vel, srate);
    // Set up uniform unity density field
    auto& dens = sim().repo().get_field("density");
    dens.setVal(rho0);
    // Set up temperature field with constant gradient in z
    auto& temp = sim().repo().get_field("temperature");
    temp.setVal(300.0);

    // Update turbulent viscosity directly
    tmodel.update_turbulent_viscosity(
        kynema_sgf::FieldState::New, DiffusionType::Crank_Nicolson);
    const auto& muturb = sim().repo().get_field("mu_turb");

    // Check values of turbulent viscosity
    auto min_val = utils::field_min(muturb);
    auto max_val = utils::field_max(muturb);
    const amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    const amrex::Real kosovic_answer =
        rho0 * kynema_sgf::utils::powi(kosovic_Cs, 2) *
        kynema_sgf::utils::powi(std::cbrt(m_dx * m_dy * m_dz), 2) * srate;
    EXPECT_NEAR(min_val, kosovic_answer, tol);
    EXPECT_NEAR(max_val, kosovic_answer, tol);

    // Check values of effective viscosity
    auto& mueff = sim().repo().get_field("velocity_mueff");
    tmodel.update_mueff(mueff);
    min_val = utils::field_min(mueff);
    max_val = utils::field_max(mueff);
    EXPECT_NEAR(min_val, kosovic_answer + 1.0e-5_rt, tol);
    EXPECT_NEAR(max_val, kosovic_answer + 1.0e-5_rt, tol);

    // Check that this effective viscosity is what gets to icns diffusion
    auto visc_name = pde_mgr.icns().fields().mueff.name();
    EXPECT_EQ(visc_name, "velocity_mueff");
}

} // namespace kynema_sgf_tests
