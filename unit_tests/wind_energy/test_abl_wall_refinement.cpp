#include "abl_test_utils.H"
#include "ks_test_utils/iter_tools.H"
#include "src/wind_energy/ABL.H"
#include "src/wind_energy/ABLWallFunction.H"
#include "src/utilities/tagging/CartBoxRefinement.H"

#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParReduce.H"
#include "AMReX_iMultiFab.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

//! Fill a field with offset + slope * log(z / z0) per component, a function
//! of the wall-normal coordinate only. The ghost cells below the wall repeat
//! the first cell so that every plane average reads finite values.
void init_log_profile(
    kynema_sgf::Field& fld,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> offset,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> slope,
    const amrex::Real z0)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();
    const int ncomp = fld.num_comp();
    for (int lev = 0; lev < nlevels; ++lev) {
        const amrex::Real dz = mesh.Geom(lev).CellSizeArray()[2];
        const amrex::Real zlo = mesh.Geom(lev).ProbLoArray()[2];
        const auto& farrs = fld(lev).arrays();
        amrex::ParallelFor(
            fld(lev), fld.num_grow(), ncomp,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                const int kk = amrex::max(k, 0);
                const amrex::Real z = zlo + ((kk + 0.5_rt) * dz);
                farrs[nbx](i, j, k, n) =
                    offset[n] + (slope[n] * std::log(z / z0));
            });
    }
    amrex::Gpu::streamSynchronize();
}

//! Mean, over the wall-adjacent cells of a level that are not covered by a
//! finer level, of the wall flux implied by the wall-model ghost value:
//! ghost * mueff / rho of the first cell
amrex::Real wall_flux_mean(
    const kynema_sgf::Field& fld,
    const kynema_sgf::Field& mueff,
    const kynema_sgf::Field& rho,
    const int comp,
    const int lev,
    amrex::Real& ncells)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();
    const int kwall = mesh.Geom(lev).Domain().smallEnd(2);

    amrex::iMultiFab level_mask;
    if (lev < nlevels - 1) {
        level_mask = amrex::makeFineMask(
            mesh.boxArray(lev), mesh.DistributionMap(lev),
            mesh.boxArray(lev + 1), mesh.refRatio(lev), 1, 0);
    } else {
        level_mask.define(
            mesh.boxArray(lev), mesh.DistributionMap(lev), 1, 0,
            amrex::MFInfo());
        level_mask.setVal(1);
    }

    const auto& f_arrs = fld(lev).const_arrays();
    const auto& mu_arrs = mueff(lev).const_arrays();
    const auto& rho_arrs = rho(lev).const_arrays();
    const auto& mask_arrs = level_mask.const_arrays();

    using SumTuple = amrex::GpuTuple<amrex::Real, amrex::Real>;
    const SumTuple sums = amrex::ParReduce(
        amrex::TypeList<amrex::ReduceOpSum, amrex::ReduceOpSum>{},
        amrex::TypeList<amrex::Real, amrex::Real>{}, fld(lev),
        amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) -> SumTuple {
            if (k != kwall) {
                return {0.0_rt, 0.0_rt};
            }
            const auto msk =
                static_cast<amrex::Real>(mask_arrs[box_no](i, j, k));
            const amrex::Real flux = f_arrs[box_no](i, j, k - 1, comp) *
                                     mu_arrs[box_no](i, j, k) /
                                     rho_arrs[box_no](i, j, k);
            return {msk * flux, msk};
        });
    amrex::GpuArray<amrex::Real, 2> vals{
        amrex::get<0>(sums), amrex::get<1>(sums)};
    amrex::ParallelDescriptor::ReduceRealSum(vals.data(), 2);
    ncells = vals[1];
    return vals[0] / amrex::max(vals[1], 1.0_rt);
}

} // namespace

/** ABL mesh with a second level covering part of the wall-modeled boundary
 *
 *  The domain is 120 x 120 x 1000 with 8 x 8 x 64 cells on level 0; level 1
 *  covers x < 60 and z < 250, so half of the wall is seen by level-1 cells
 *  whose centers sit at half the height of the level-0 wall cells.
 */
class ABLWallRefinementTest : public ABLMeshTest
{
protected:
    void populate_parameters() override
    {
        ABLMeshTest::populate_parameters();
        {
            amrex::ParmParse pp("amr");
            pp.add("max_level", 1);
            pp.add("max_grid_size", 64);
            pp.add("blocking_factor", 4);
            pp.add("n_error_buf", 0);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<int> periodic{{1, 1, 0}};
            pp.addarr("is_periodic", periodic);
        }
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
        {
            amrex::ParmParse pp("incflo");
            pp.add("diffusion_type", 0);
        }
        {
            amrex::ParmParse pp("transport");
            pp.add("viscosity", m_mu);
            pp.add("laminar_prandtl", 1.0_rt);
        }
        {
            amrex::ParmParse pp("ABL");
            pp.add("kappa", m_kappa);
            pp.add("surface_roughness_z0", m_z0);
            pp.add("surface_temp_flux", m_qwall);
            pp.add("wall_shear_stress_type", m_shear_stress_type);
        }

        // Level-1 box over x < 60, all y, z < 250
        std::stringstream ss;
        ss << "1 // Number of levels" << '\n';
        ss << "1 // Number of boxes at this level" << '\n';
        ss << "0.0 0.0 0.0 60.0 120.0 250.0" << '\n';

        create_mesh_instance<RefineMesh>();
        auto box_refine =
            std::make_unique<kynema_sgf::CartBoxRefinement>(sim());
        box_refine->read_inputs(mesh(), ss);
        if (mesh<RefineMesh>() != nullptr) {
            mesh<RefineMesh>()->refine_criteria_vec().push_back(
                std::move(box_refine));
        }
    }

    //! Check that the mean wall stress and heat flux of every level match
    //! the friction velocity and the surface heat flux of the wall function
    void check_level_fluxes()
    {
        constexpr amrex::Real tol =
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
        constexpr amrex::Real ustar = 0.5_rt;
        constexpr amrex::Real thetastar = -0.05_rt;
        constexpr amrex::Real wind_angle = 0.5_rt;
        constexpr amrex::Real theta0 = 300.0_rt;

        populate_parameters();
        initialize_mesh();
        ASSERT_EQ(sim().repo().num_active_levels(), 2);

        auto& pde_mgr = sim().pde_manager();
        pde_mgr.register_icns();
        sim().create_turbulence_model();
        sim().init_physics();

        auto& repo = sim().repo();
        auto& velocity = repo.get_field("velocity");
        auto& temperature = repo.get_field("temperature");
        auto& density = repo.get_field("density");
        auto& vel_mueff = repo.get_field("velocity_mueff");
        auto& temp_mueff = repo.get_field("temperature_mueff");

        // Logarithmic wind and temperature profiles, uniform in every plane
        const amrex::Real wind_cos = std::cos(wind_angle);
        const amrex::Real wind_sin = std::sin(wind_angle);
        init_log_profile(
            velocity, {0.0_rt, 0.0_rt, 0.0_rt},
            {ustar / m_kappa * wind_cos, ustar / m_kappa * wind_sin, 0.0_rt},
            m_z0);
        init_log_profile(
            temperature, {theta0, 0.0_rt, 0.0_rt},
            {thetastar / m_kappa, 0.0_rt, 0.0_rt}, m_z0);
        density.setVal(1.0_rt);
        vel_mueff.setVal(m_mu);
        temp_mueff.setVal(m_mu);

        // Wall function: plane averages, friction velocity, custom BCs
        for (auto& pp : sim().physics()) {
            pp->post_init_actions();
        }
        pde_mgr.advance_states();

        // Fill the wall-model ghost cells of velocity and temperature
        velocity.apply_bc_funcs(kynema_sgf::FieldState::Old);
        temperature.apply_bc_funcs(kynema_sgf::FieldState::Old);

        const auto& abl = sim().physics_manager().get<kynema_sgf::ABL>();
        const auto& mo = abl.abl_wall_function().mo();
        ASSERT_GT(mo.utau, 0.0_rt);
        const amrex::Real utau2 = mo.utau * mo.utau;
        // The heat flux is a small difference of temperatures of order
        // theta0, so its roundoff scales with theta0
        const amrex::Real tol_q = tol * theta0 * mo.utau;

        // Every level must return the same mean stress, u_*^2 along the mean
        // wind, and the same mean heat flux, the specified surface flux (the
        // ghost value is the wall-normal gradient, so the flux is minus it)
        for (int lev = 0; lev < 2; ++lev) {
            amrex::Real ncells = 0.0_rt;
            const amrex::Real taux =
                wall_flux_mean(velocity, vel_mueff, density, 0, lev, ncells);
            // Both levels must own part of the wall
            ASSERT_GT(ncells, 0.0_rt) << "level " << lev;
            const amrex::Real tauy =
                wall_flux_mean(velocity, vel_mueff, density, 1, lev, ncells);
            const amrex::Real qwall = -wall_flux_mean(
                temperature, temp_mueff, density, 0, lev, ncells);
            EXPECT_NEAR(taux, utau2 * wind_cos, tol * utau2) << "level " << lev;
            EXPECT_NEAR(tauy, utau2 * wind_sin, tol * utau2) << "level " << lev;
            EXPECT_NEAR(qwall, m_qwall, tol_q) << "level " << lev;
        }
    }

    std::string m_shear_stress_type{"moeng"};
    const amrex::Real m_mu{0.01_rt};
    const amrex::Real m_kappa{0.41_rt};
    const amrex::Real m_z0{0.1_rt};
    const amrex::Real m_qwall{0.02_rt};
};

TEST_F(ABLWallRefinementTest, moeng_wall_model_level_consistent)
{
    m_shear_stress_type = "moeng";
    check_level_fluxes();
}

TEST_F(ABLWallRefinementTest, schumann_wall_model_level_consistent)
{
    m_shear_stress_type = "schumann";
    check_level_fluxes();
}

} // namespace kynema_sgf_tests
