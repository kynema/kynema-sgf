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

//! Multiply the valid cells of a level that are covered by the next finer
//! level by a factor
void scale_covered_cells(
    kynema_sgf::Field& fld, const int lev, const amrex::Real factor)
{
    const auto& mesh = fld.repo().mesh();
    const auto mask = amrex::makeFineMask(
        fld(lev), mesh.boxArray(lev + 1), mesh.refRatio(lev), 0, 1);
    const auto& farrs = fld(lev).arrays();
    const auto& mask_arrs = mask.const_arrays();
    amrex::ParallelFor(
        fld(lev), amrex::IntVect(0), fld.num_comp(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
            if (mask_arrs[nbx](i, j, k) == 1) {
                farrs[nbx](i, j, k, n) *= factor;
            }
        });
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

//! Add a perturbation that varies across the wall plane, so that a mean
//! taken over the wrong cells or boxes cannot come out right by symmetry
void perturb_in_plane(kynema_sgf::Field& fld, const amrex::Real amp)
{
    const auto& mesh = fld.repo().mesh();
    const int nlevels = fld.repo().num_active_levels();
    const int ncomp = fld.num_comp();
    for (int lev = 0; lev < nlevels; ++lev) {
        const int rr = (lev == 0) ? 1 : mesh.refRatio(0)[0];
        const auto& farrs = fld(lev).arrays();
        amrex::ParallelFor(
            fld(lev), fld.num_grow(), ncomp,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                // In level-0 index units, so levels see the same pattern
                farrs[nbx](i, j, k, n) +=
                    amp * (1.0_rt + n) *
                    (static_cast<amrex::Real>(i / rr) +
                     (7.0_rt * static_cast<amrex::Real>(j / rr)));
            });
    }
    amrex::Gpu::streamSynchronize();
}

//! Blank the wall cells of a level with i below a bound, as terrain would,
//! and set the velocity there to rest
void blank_wall_cells(
    amrex::iMultiFab& blank,
    kynema_sgf::Field& vel,
    const int kwall,
    const int ibound)
{
    const auto& blank_arrs = blank.arrays();
    const auto& vel_arrs = vel(1).arrays();
    amrex::ParallelFor(
        blank, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const int b = ((k == kwall) && (i < ibound)) ? 1 : 0;
            blank_arrs[nbx](i, j, k) = b;
            if (b == 1) {
                vel_arrs[nbx](i, j, k, 0) = 0.0_rt;
                vel_arrs[nbx](i, j, k, 1) = 0.0_rt;
            }
        });
    amrex::Gpu::streamSynchronize();
}

/** Mean quantities of the uncovered wall cells of a level, summed over every
 *  cell of the level (the reference the wall-layer reduction must match)
 *
 *  \return u, v, |u_h|, |u_h| u, |u_h| v and theta means
 */
amrex::GpuArray<amrex::Real, 6> full_level_wall_means(
    const kynema_sgf::Field& vel, const kynema_sgf::Field& theta, const int lev)
{
    const auto& mesh = vel.repo().mesh();
    const int nlevels = vel.repo().num_active_levels();
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
    const auto& v_arrs = vel(lev).const_arrays();
    const auto& t_arrs = theta(lev).const_arrays();
    const auto& mask_arrs = level_mask.const_arrays();
    using Tuple = amrex::GpuTuple<
        amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real,
        amrex::Real, amrex::Real>;
    const Tuple sums = amrex::ParReduce(
        amrex::TypeList<
            amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
            amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
            amrex::ReduceOpSum>{},
        amrex::TypeList<
            amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real,
            amrex::Real, amrex::Real>{},
        vel(lev), amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) -> Tuple {
            if (k != kwall) {
                return {0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt, 0.0_rt};
            }
            const auto m = static_cast<amrex::Real>(mask_arrs[box_no](i, j, k));
            const amrex::Real uu = v_arrs[box_no](i, j, k, 0);
            const amrex::Real vv = v_arrs[box_no](i, j, k, 1);
            const amrex::Real ws = std::sqrt((uu * uu) + (vv * vv));
            return {m * uu,
                    m * vv,
                    m * ws,
                    m * ws * uu,
                    m * ws * vv,
                    m * t_arrs[box_no](i, j, k),
                    m};
        });
    amrex::GpuArray<amrex::Real, 7> v{amrex::get<0>(sums), amrex::get<1>(sums),
                                      amrex::get<2>(sums), amrex::get<3>(sums),
                                      amrex::get<4>(sums), amrex::get<5>(sums),
                                      amrex::get<6>(sums)};
    amrex::ParallelDescriptor::ReduceRealSum(v.data(), 7);
    const amrex::Real n = amrex::max(v[6], 1.0_rt);
    return {v[0] / n, v[1] / n, v[2] / n, v[3] / n, v[4] / n, v[5] / n};
}

} // namespace

/** A RefineMesh whose levels list their boxes in reverse order
 *
 *  AMReX lists the boxes on the wall first, so the boxes of the wall layer
 *  of a level normally have the same indices as the level's own boxes.
 *  Reversing the order makes them differ, so a test can tell whether the
 *  wall layer reads the right box of each field.
 */
class ReversedRefineMesh : public RefineMesh
{
protected:
    void MakeNewLevelFromScratch(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& /*dm*/) override
    {
        amrex::BoxList bl;
        for (int ib = static_cast<int>(ba.size()) - 1; ib >= 0; --ib) {
            bl.push_back(ba[ib]);
        }
        const amrex::BoxArray rba(std::move(bl));
        const amrex::DistributionMapping rdm(rba);
        RefineMesh::MakeNewLevelFromScratch(lev, time, rba, rdm);
    }
};

/** ABL mesh with a second level covering part of the wall-modeled boundary
 *
 *  The domain is 120 x 120 x 1000 with 8 x 8 x 64 cells on level 0; level 1
 *  is tagged over x < 60 and z < 250 and, rounded out to the blocking
 *  factor, covers x < 90 and z < 281.25, so part of the wall is seen by
 *  level-1 cells
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
            pp.add("max_grid_size", m_max_grid_size);
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
            if (m_wall_position_set) {
                pp.add("wall_position", m_wall_position);
            }
            pp.add("surface_temp_flux", m_qwall);
            pp.add("wall_shear_stress_type", m_shear_stress_type);
            if (m_log_law_height > 0.0_rt) {
                pp.add("log_law_height", m_log_law_height);
            }
        }

        std::stringstream ss;
        ss << "1 // Number of levels" << '\n';
        ss << "1 // Number of boxes at this level" << '\n';
        ss << m_refine_box << '\n';

        if (m_reverse_boxes) {
            create_mesh_instance<ReversedRefineMesh>();
        } else {
            create_mesh_instance<RefineMesh>();
        }
        auto box_refine =
            std::make_unique<kynema_sgf::CartBoxRefinement>(sim());
        box_refine->read_inputs(mesh(), ss);
        if (mesh<RefineMesh>() != nullptr) {
            mesh<RefineMesh>()->refine_criteria_vec().push_back(
                std::move(box_refine));
        }
    }

    //! Set up the two-level mesh with plane-uniform logarithmic profiles and
    //! fill the wall-model ghost cells of velocity and temperature
    void init_wall_fields()
    {
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
        const amrex::Real wind_cos = std::cos(m_wind_angle);
        const amrex::Real wind_sin = std::sin(m_wind_angle);
        init_log_profile(
            velocity, {0.0_rt, 0.0_rt, 0.0_rt},
            {m_ustar / m_kappa * wind_cos, m_ustar / m_kappa * wind_sin,
             0.0_rt},
            m_z0);
        init_log_profile(
            temperature, {m_theta0, 0.0_rt, 0.0_rt},
            {m_thetastar / m_kappa, 0.0_rt, 0.0_rt}, m_z0);
        if (m_blank_ibound > 0) {
            // Terrain blanking half of the level-1 wall cells (level 1
            // covers i < 16 there)
            auto& blank = repo.declare_int_field("terrain_blank", 1, 1, 1);
            blank.setVal(0);
            blank_wall_cells(
                blank(1), velocity,
                sim().repo().mesh().Geom(1).Domain().smallEnd(2),
                m_blank_ibound);
        }
        if (m_perturb != 0.0_rt) {
            perturb_in_plane(velocity, m_perturb);
            perturb_in_plane(temperature, m_perturb);
        }
        if (m_covered_scale != 1.0_rt) {
            // Level-0 cells covered by level 1, which the means of level 0
            // must not see
            scale_covered_cells(velocity, 0, m_covered_scale);
            scale_covered_cells(temperature, 0, m_covered_scale);
        }
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
    }

    //! Check that the mean wall stress and heat flux of every level match
    //! the friction velocity and the surface heat flux of the wall function
    void check_level_fluxes()
    {
        constexpr amrex::Real tol =
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

        ASSERT_NO_FATAL_FAILURE(init_wall_fields());

        auto& repo = sim().repo();
        const auto& velocity = repo.get_field("velocity");
        const auto& temperature = repo.get_field("temperature");
        const auto& density = repo.get_field("density");
        const auto& vel_mueff = repo.get_field("velocity_mueff");
        const auto& temp_mueff = repo.get_field("temperature_mueff");
        const amrex::Real wind_cos = std::cos(m_wind_angle);
        const amrex::Real wind_sin = std::sin(m_wind_angle);

        const auto& abl = sim().physics_manager().get<kynema_sgf::ABL>();
        const auto& mo = abl.abl_wall_function().mo();
        ASSERT_GT(mo.utau, 0.0_rt);
        const amrex::Real utau2 = mo.utau * mo.utau;
        // The heat flux is a small difference of temperatures of order
        // theta0, so its roundoff scales with theta0
        const amrex::Real tol_q = 1.0e2_rt *
                                  std::numeric_limits<amrex::Real>::epsilon() *
                                  m_theta0 * mo.utau;

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

    //! Check that the Donelan model uses the mean wind at the reference
    //! height on every level: the stress of the wall cells of each level is
    //! Cd(U_ref) |u_h| u_h with the drag coefficient of the reference-height
    //! mean wind, not of the first-cell mean wind of that level
    void check_donelan_reference_height()
    {
        constexpr amrex::Real tol =
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

        ASSERT_NO_FATAL_FAILURE(init_wall_fields());

        auto& repo = sim().repo();
        const auto& velocity = repo.get_field("velocity");
        const auto& density = repo.get_field("density");
        const auto& vel_mueff = repo.get_field("velocity_mueff");
        const amrex::Real wind_cos = std::cos(m_wind_angle);
        const amrex::Real wind_sin = std::sin(m_wind_angle);

        const auto& abl = sim().physics_manager().get<kynema_sgf::ABL>();
        const auto& wall_func = abl.abl_wall_function();
        // The mesh has per-level mean quantities, which Donelan must not use
        ASSERT_NE(&wall_func.mo(1), &wall_func.mo());
        const amrex::Real wspd_ref = wall_func.mo().vmag_mean;
        // Drag coefficient of ShearStressDonelan, in its linear range for
        // the mean wind of the reference height and of both first cells
        ASSERT_GT(wspd_ref, 5.0_rt);
        ASSERT_LT(wspd_ref, 25.0_rt);
        const amrex::Real cd = 0.001_rt + (7.0e-5_rt * (wspd_ref - 5.0_rt));

        for (int lev = 0; lev < 2; ++lev) {
            // Wind speed of the first cell of this level, uniform in plane
            const amrex::Real z1 =
                0.5_rt * repo.mesh().Geom(lev).CellSizeArray()[2];
            const amrex::Real wspd = m_ustar / m_kappa * std::log(z1 / m_z0);
            ASSERT_GT(wspd, 5.0_rt) << "level " << lev;
            const amrex::Real tau = cd * wspd * wspd;
            amrex::Real ncells = 0.0_rt;
            const amrex::Real taux =
                wall_flux_mean(velocity, vel_mueff, density, 0, lev, ncells);
            ASSERT_GT(ncells, 0.0_rt) << "level " << lev;
            const amrex::Real tauy =
                wall_flux_mean(velocity, vel_mueff, density, 1, lev, ncells);
            EXPECT_NEAR(taux, tau * wind_cos, tol * tau) << "level " << lev;
            EXPECT_NEAR(tauy, tau * wind_sin, tol * tau) << "level " << lev;
        }
    }

    //! Check that the mean quantities of every level are those of the
    //! wall-adjacent cells that the level owns, with the level-0 cells
    //! covered by level 1 set to other values
    void check_covered_cells_excluded()
    {
        constexpr amrex::Real tol =
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;

        m_covered_scale = 2.0_rt;
        ASSERT_NO_FATAL_FAILURE(init_wall_fields());

        const auto& mesh = sim().repo().mesh();
        const auto& abl = sim().physics_manager().get<kynema_sgf::ABL>();
        const auto& wall_func = abl.abl_wall_function();
        const amrex::Real wind_cos = std::cos(m_wind_angle);
        const amrex::Real wind_sin = std::sin(m_wind_angle);
        for (int lev = 0; lev < 2; ++lev) {
            const auto& mo_lev = wall_func.mo(lev);
            ASSERT_NE(&mo_lev, &wall_func.mo()) << "level " << lev;
            // Plane-uniform first-cell values of the cells the level owns
            const amrex::Real z1 = 0.5_rt * mesh.Geom(lev).CellSizeArray()[2];
            const amrex::Real wspd = m_ustar / m_kappa * std::log(z1 / m_z0);
            const amrex::Real theta =
                m_theta0 + (m_thetastar / m_kappa * std::log(z1 / m_z0));
            EXPECT_NEAR(mo_lev.zref, z1, tol * z1) << "level " << lev;
            EXPECT_NEAR(mo_lev.vel_mean[0], wspd * wind_cos, tol * wspd)
                << "level " << lev;
            EXPECT_NEAR(mo_lev.vel_mean[1], wspd * wind_sin, tol * wspd)
                << "level " << lev;
            EXPECT_NEAR(mo_lev.vmag_mean, wspd, tol * wspd) << "level " << lev;
            EXPECT_NEAR(
                mo_lev.Su_mean, wspd * wspd * wind_cos, tol * wspd * wspd)
                << "level " << lev;
            EXPECT_NEAR(
                mo_lev.Sv_mean, wspd * wspd * wind_sin, tol * wspd * wspd)
                << "level " << lev;
            EXPECT_NEAR(
                mo_lev.theta_mean, theta,
                1.0e2_rt * std::numeric_limits<amrex::Real>::epsilon() *
                    m_theta0)
                << "level " << lev;
        }
    }

    //! Level-1 box (xlo ylo zlo xhi yhi zhi), by default over x < 60, all
    //! y and z < 250
    std::string m_refine_box{"0.0 0.0 0.0 60.0 120.0 250.0"};
    std::string m_shear_stress_type{"moeng"};
    //! Factor on the level-0 cells covered by level 1 (1: plane uniform)
    amrex::Real m_covered_scale{1.0_rt};
    //! Amplitude of an in-plane perturbation of the fields (0: none)
    amrex::Real m_perturb{0.0_rt};
    //! Level-1 wall cells with i below this are blanked by terrain (0: none)
    int m_blank_ibound{0};
    //! ABL.wall_position, when set
    bool m_wall_position_set{false};
    amrex::Real m_wall_position{0.0_rt};
    //! Largest box size; small values give several boxes per level
    int m_max_grid_size{64};
    //! List the boxes of every level in reverse order
    bool m_reverse_boxes{false};
    //! Reference height of the plane averages, first-cell height if <= 0
    amrex::Real m_log_law_height{0.0_rt};
    amrex::Real m_ustar{0.5_rt};
    const amrex::Real m_thetastar{-0.05_rt};
    const amrex::Real m_wind_angle{0.5_rt};
    const amrex::Real m_theta0{300.0_rt};
    const amrex::Real m_mu{0.01_rt};
    const amrex::Real m_kappa{0.41_rt};
    amrex::Real m_z0{0.1_rt};
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

TEST_F(ABLWallRefinementTest, local_wall_model_level_consistent)
{
    m_shear_stress_type = "local";
    check_level_fluxes();
}

// The constant model uses the plane averages only; this guards that each
// level still returns the friction velocity and the surface heat flux
TEST_F(ABLWallRefinementTest, constant_wall_model_level_consistent)
{
    m_shear_stress_type = "constant";
    check_level_fluxes();
}

TEST_F(ABLWallRefinementTest, donelan_wall_model_keeps_reference_height)
{
    // Donelan selects its drag coefficient from the mean wind at the
    // reference height, here 10 m as in the hurricane boundary layer case,
    // with a wind strong enough for the linear range of the coefficient
    m_shear_stress_type = "donelan";
    m_log_law_height = 10.0_rt;
    m_ustar = 1.0_rt;
    check_donelan_reference_height();
}

TEST_F(ABLWallRefinementTest, level_means_exclude_covered_wall_cells)
{
    check_covered_cells_excluded();
}

// With several boxes per level, listed so that the wall boxes are not the
// first ones (and, on several ranks, ranks that own no wall box), and fields
// that vary across the plane, the level means of the wall-layer reduction
// equal a sum over every cell of the level
TEST_F(ABLWallRefinementTest, level_means_match_a_full_sum_on_many_boxes)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    m_max_grid_size = 4;
    m_perturb = 1.0e-2_rt;
    m_reverse_boxes = true;
    ASSERT_NO_FATAL_FAILURE(init_wall_fields());
    // The wall boxes must not be the first ones, or the wall layer would
    // have the same box indices as the level whatever the code does
    const auto& ba0 = sim().repo().mesh().boxArray(0);
    ASSERT_GT(ba0.size(), 4);
    ASSERT_GT(
        ba0[0].smallEnd(2), sim().repo().mesh().Geom(0).Domain().smallEnd(2));

    auto& repo = sim().repo();
    const auto& abl = sim().physics_manager().get<kynema_sgf::ABL>();
    const auto& wall_func = abl.abl_wall_function();
    for (int lev = 0; lev < 2; ++lev) {
        const auto& mo_lev = wall_func.mo(lev);
        ASSERT_NE(&mo_lev, &wall_func.mo()) << "level " << lev;
        const auto ref = full_level_wall_means(
            repo.get_field("velocity"), repo.get_field("temperature"), lev);
        EXPECT_NEAR(mo_lev.vel_mean[0], ref[0], tol * std::abs(ref[2]))
            << "level " << lev;
        EXPECT_NEAR(mo_lev.vel_mean[1], ref[1], tol * std::abs(ref[2]))
            << "level " << lev;
        EXPECT_NEAR(mo_lev.vmag_mean, ref[2], tol * std::abs(ref[2]))
            << "level " << lev;
        EXPECT_NEAR(mo_lev.Su_mean, ref[3], tol * ref[2] * ref[2])
            << "level " << lev;
        EXPECT_NEAR(mo_lev.Sv_mean, ref[4], tol * ref[2] * ref[2])
            << "level " << lev;
        EXPECT_NEAR(mo_lev.theta_mean, ref[5], tol * m_theta0)
            << "level " << lev;
    }
}

// Wall cells blanked by terrain carry no wall stress, so the level means are
// those of the fluid cells: here half the level-1 wall cells are blanked and
// at rest, and the level-1 mean wind is still the profile's first-cell wind
TEST_F(ABLWallRefinementTest, level_means_exclude_blanked_wall_cells)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    m_blank_ibound = 8;
    ASSERT_NO_FATAL_FAILURE(init_wall_fields());

    const auto& mesh = sim().repo().mesh();
    const auto& abl = sim().physics_manager().get<kynema_sgf::ABL>();
    const auto& mo1 = abl.abl_wall_function().mo(1);
    ASSERT_NE(&mo1, &abl.abl_wall_function().mo());
    const amrex::Real z1 = 0.5_rt * mesh.Geom(1).CellSizeArray()[2];
    const amrex::Real wspd = m_ustar / m_kappa * std::log(z1 / m_z0);
    EXPECT_NEAR(mo1.vmag_mean, wspd, tol * wspd);
    EXPECT_NEAR(mo1.vel_mean[0], wspd * std::cos(m_wind_angle), tol * wspd);
}

// A level whose first cell sits at or below the roughness height keeps the
// means at the reference height: the log law, and phi_h, do not hold there
TEST_F(ABLWallRefinementTest, first_cell_below_roughness_keeps_reference)
{
    // Level-1 first cells at 3.9 m, level-0 ones at 7.8 m
    m_z0 = 5.0_rt;
    ASSERT_NO_FATAL_FAILURE(init_wall_fields());

    const auto& abl = sim().physics_manager().get<kynema_sgf::ABL>();
    const auto& wall_func = abl.abl_wall_function();
    EXPECT_EQ(wall_func.mo(1).zref, wall_func.mo().zref);
    EXPECT_EQ(wall_func.mo(1).vmag_mean, wall_func.mo().vmag_mean);
    EXPECT_GT(wall_func.mo(0).zref, m_z0);
}

// The first-cell height of a level is measured from ABL.wall_position, like
// the reference height
TEST_F(ABLWallRefinementTest, level_height_is_measured_from_the_wall)
{
    constexpr amrex::Real tol =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt;
    m_wall_position_set = true;
    m_wall_position = -2.0_rt;
    ASSERT_NO_FATAL_FAILURE(init_wall_fields());

    const auto& mesh = sim().repo().mesh();
    const auto& abl = sim().physics_manager().get<kynema_sgf::ABL>();
    const auto& mo1 = abl.abl_wall_function().mo(1);
    ASSERT_NE(&mo1, &abl.abl_wall_function().mo());
    const amrex::Real z1 =
        (0.5_rt * mesh.Geom(1).CellSizeArray()[2]) - m_wall_position;
    EXPECT_NEAR(mo1.zref, z1, tol * z1);
}

TEST_F(ABLWallRefinementTest, refinement_aloft_keeps_reference_height)
{
    // Level-1 box that does not reach the wall: level 0 alone owns the
    // wall and keeps the plane averages at the reference height, as on a
    // single-level mesh
    m_refine_box = "0.0 0.0 500.0 60.0 120.0 750.0";
    ASSERT_NO_FATAL_FAILURE(init_wall_fields());

    const auto& abl = sim().physics_manager().get<kynema_sgf::ABL>();
    const auto& wall_func = abl.abl_wall_function();
    EXPECT_EQ(&wall_func.mo(0), &wall_func.mo());
    EXPECT_EQ(&wall_func.mo(1), &wall_func.mo());
}

} // namespace kynema_sgf_tests
