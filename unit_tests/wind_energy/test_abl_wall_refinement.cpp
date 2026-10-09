#include "abl_test_utils.H"
#include "src/wind_energy/ABL.H"
#include "src/wind_energy/ABLWallFunction.H"
#include "src/wind_energy/MOData.H"
#include "src/utilities/trig_ops.H"
#include "src/utilities/tagging/CartBoxRefinement.H"

#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParReduce.H"
#include "AMReX_iMultiFab.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {

//! Monin-Obukhov surface-layer profile (Dyer functions, as in MOData) for a
//! friction velocity, a temperature scale and a surface temperature. The
//! Obukhov length is the fixed point of MOData::update_fluxes fed with the
//! exact plane means at the reference height, so that the wall function
//! recovers u*, L and q = -u* theta* from the profile (to its own convergence
//! tolerance of 1e-5 on u*).
struct MOProfile
{
    amrex::Real ustar{0.5_rt};
    amrex::Real thetastar{-0.04_rt};
    amrex::Real theta_s{300.0_rt};
    amrex::Real z0{0.1_rt};
    amrex::Real kappa{0.41_rt};
    amrex::Real obl{std::numeric_limits<amrex::Real>::max()};

    //! Upward surface heat flux
    [[nodiscard]] amrex::Real q() const { return -ustar * thetastar; }

    [[nodiscard]] amrex::Real phi_m(const amrex::Real z) const
    {
        return std::log(z / z0) -
               kynema_sgf::MOData::calc_psi_m(z / obl, 16.0_rt, 5.0_rt);
    }

    [[nodiscard]] amrex::Real phi_h(const amrex::Real z) const
    {
        return std::log(z / z0) -
               kynema_sgf::MOData::calc_psi_h(z / obl, 16.0_rt, 5.0_rt);
    }

    [[nodiscard]] amrex::Real speed(const amrex::Real z) const
    {
        return ustar / kappa * phi_m(z);
    }

    [[nodiscard]] amrex::Real theta(const amrex::Real z) const
    {
        return theta_s + (thetastar / kappa * phi_h(z));
    }

    //! Obukhov length consistent with update_fluxes at the reference height
    void solve_obukhov_length(const amrex::Real zref)
    {
        obl = std::numeric_limits<amrex::Real>::max();
        if (std::abs(q()) < std::numeric_limits<amrex::Real>::epsilon()) {
            return;
        }
        for (int it = 0; it < 50; ++it) {
            obl =
                -ustar * ustar * ustar * theta(zref) / (kappa * 9.81_rt * q());
        }
    }
};

//! Fill velocity and temperature with the Monin-Obukhov profile, uniform in
//! every plane, with heights measured from the wall position. The ghost cells
//! below the wall repeat the first cell so that every plane average reads
//! finite values
void init_mo_profiles(
    kynema_sgf::Field& vel,
    kynema_sgf::Field& temp,
    const MOProfile& prof,
    const amrex::Real wind_angle,
    const amrex::Real wall_pos)
{
    const amrex::Real ustar = prof.ustar;
    const amrex::Real thetastar = prof.thetastar;
    const amrex::Real theta_s = prof.theta_s;
    const amrex::Real z0 = prof.z0;
    const amrex::Real kappa = prof.kappa;
    const amrex::Real obl = prof.obl;
    const amrex::Real wc = std::cos(wind_angle);
    const amrex::Real ws = std::sin(wind_angle);
    const auto& mesh = vel.repo().mesh();
    const int nlevels = vel.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const amrex::Real dz = mesh.Geom(lev).CellSizeArray()[2];
        const amrex::Real zlo = mesh.Geom(lev).ProbLoArray()[2];
        const auto& varrs = vel(lev).arrays();
        const auto& tarrs = temp(lev).arrays();
        const amrex::IntVect ng = amrex::min(vel.num_grow(), temp.num_grow());
        amrex::ParallelFor(
            vel(lev), ng, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const int kk = amrex::max(k, 0);
                const amrex::Real z = zlo + ((kk + 0.5_rt) * dz) - wall_pos;
                const amrex::Real zeta = z / obl;
                amrex::Real psim = -5.0_rt * zeta;
                amrex::Real psih = -5.0_rt * zeta;
                if (zeta < 0.0_rt) {
                    const amrex::Real xm =
                        std::sqrt(std::sqrt(1.0_rt - (16.0_rt * zeta)));
                    psim = (2.0_rt * std::log(0.5_rt * (1.0_rt + xm))) +
                           std::log(0.5_rt * (1.0_rt + (xm * xm))) -
                           (2.0_rt * std::atan(xm)) +
                           kynema_sgf::utils::half_pi();
                    const amrex::Real xh = std::sqrt(1.0_rt - (16.0_rt * zeta));
                    psih = 2.0_rt * std::log(0.5_rt * (1.0_rt + xh));
                }
                const amrex::Real spd =
                    ustar / kappa * (std::log(z / z0) - psim);
                varrs[nbx](i, j, k, 0) = spd * wc;
                varrs[nbx](i, j, k, 1) = spd * ws;
                varrs[nbx](i, j, k, 2) = 0.0_rt;
                tarrs[nbx](i, j, k, 0) =
                    theta_s + (thetastar / kappa * (std::log(z / z0) - psih));
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

/** ABL mesh with a second level covering part or all of the wall-modeled
 *  boundary
 *
 *  The domain is 120 x 120 x 1000 with 8 x 8 x 64 cells on level 0 (first
 *  cell at 7.8125 m); level 1 (first cell at 3.90625 m) is tagged over
 *  x < 60 and z < 250 and, rounded out to the blocking factor, covers x < 90
 *  and z < 281.25, so three quarters of the wall are seen by level-1 cells.
 *  With m_full_coverage the level-1 box covers the whole wall.
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
            pp.add("kappa", m_prof.kappa);
            pp.add("surface_roughness_z0", m_prof.z0);
            if (m_wall_position_set) {
                pp.add("wall_position", m_wall_position);
            }
            if (m_surface_temp_mode) {
                // Surface-temperature mode: a fixed surface temperature
                // equal to the profile's
                pp.add("surface_temp_rate", 0.0_rt);
                pp.add("surface_temp_init", m_prof.theta_s);
            } else {
                pp.add("surface_temp_flux", m_prof.q());
            }
            pp.add("wall_shear_stress_type", m_shear_stress_type);
            if (m_log_law_height > 0.0_rt) {
                pp.add("log_law_height", m_log_law_height);
            }
            if (!m_level_means_method.empty()) {
                pp.add("level_means_method", m_level_means_method);
            }
        }

        std::stringstream ss;
        ss << "1 // Number of levels" << '\n';
        ss << "1 // Number of boxes at this level" << '\n';
        ss << (m_full_coverage ? "0.0 0.0 0.0 120.0 120.0 250.0" : m_refine_box)
           << '\n';

        create_mesh_instance<RefineMesh>();
        auto box_refine =
            std::make_unique<kynema_sgf::CartBoxRefinement>(sim());
        box_refine->read_inputs(mesh(), ss);
        mesh<RefineMesh>()->refine_criteria_vec().push_back(
            std::move(box_refine));
    }

    //! Mesh, equation systems and physics (the wall function reads its
    //! inputs here)
    void init_physics()
    {
        populate_parameters();
        initialize_mesh();
        ASSERT_EQ(sim().repo().num_active_levels(), 2);
        sim().pde_manager().register_icns();
        sim().create_turbulence_model();
        sim().init_physics();
    }

    //! Set up the two-level mesh with plane-uniform Monin-Obukhov profiles
    //! and fill the wall-model ghost cells of velocity and temperature
    void init_wall_fields()
    {
        ASSERT_NO_FATAL_FAILURE(init_physics());

        auto& repo = sim().repo();
        auto& velocity = repo.get_field("velocity");
        auto& temperature = repo.get_field("temperature");
        auto& density = repo.get_field("density");
        auto& vel_mueff = repo.get_field("velocity_mueff");
        auto& temp_mueff = repo.get_field("temperature_mueff");

        // The profile's Obukhov length is the one the wall function finds
        // from the exact plane means at the reference height
        const amrex::Real zref =
            (m_log_law_height > 0.0_rt)
                ? m_log_law_height
                : 0.5_rt * repo.mesh().Geom(0).CellSizeArray()[2];
        m_prof.solve_obukhov_length(zref);
        const amrex::Real wall_pos =
            m_wall_position_set ? m_wall_position : 0.0_rt;
        init_mo_profiles(velocity, temperature, m_prof, m_wind_angle, wall_pos);
        density.setVal(1.0_rt);
        vel_mueff.setVal(m_mu);
        temp_mueff.setVal(m_mu);

        // Wall function: plane averages, friction velocity, custom BCs
        for (auto& pp : sim().physics()) {
            pp->post_init_actions();
        }
        sim().pde_manager().advance_states();

        // Fill the wall-model ghost cells of velocity and temperature
        velocity.apply_bc_funcs(kynema_sgf::FieldState::Old);
        temperature.apply_bc_funcs(kynema_sgf::FieldState::Old);
    }

    //! Mean wall stress along the mean wind, its crosswind part and the mean
    //! upward heat flux of the wall cells a level owns (the ghost value is
    //! the wall-normal gradient, so the flux is minus it)
    void level_fluxes(
        const int lev,
        amrex::Real& tau_along,
        amrex::Real& tau_cross,
        amrex::Real& qwall,
        amrex::Real& ncells)
    {
        auto& repo = sim().repo();
        const auto& velocity = repo.get_field("velocity");
        const auto& temperature = repo.get_field("temperature");
        const auto& density = repo.get_field("density");
        const auto& vel_mueff = repo.get_field("velocity_mueff");
        const auto& temp_mueff = repo.get_field("temperature_mueff");
        const amrex::Real wc = std::cos(m_wind_angle);
        const amrex::Real ws = std::sin(m_wind_angle);
        const amrex::Real taux =
            wall_flux_mean(velocity, vel_mueff, density, 0, lev, ncells);
        const amrex::Real tauy =
            wall_flux_mean(velocity, vel_mueff, density, 1, lev, ncells);
        tau_along = (taux * wc) + (tauy * ws);
        tau_cross = (tauy * wc) - (taux * ws);
        qwall =
            -wall_flux_mean(temperature, temp_mueff, density, 0, lev, ncells);
    }

    //! Check that every level with wall cells returns the profile's friction
    //! velocity squared along the wind and the profile's heat flux. The
    //! reference height is set at the centre of a level-1 cell and level 1
    //! covers the whole wall, so the plane averages the wall function reads
    //! are the exact profile values. The friction-velocity iteration of
    //! MOData::update_fluxes converges to about 1e-7; the tolerance leaves a
    //! margin for single precision
    void check_level_fluxes_exact_reference()
    {
        constexpr amrex::Real tol = 1.0e-4_rt;
        m_full_coverage = true;
        m_log_law_height = 11.71875_rt; // centre of the second level-1 cell
        ASSERT_NO_FATAL_FAILURE(init_wall_fields());

        const amrex::Real utau2 = m_prof.ustar * m_prof.ustar;
        // The heat flux is a small difference of temperatures of order
        // theta_s, so its roundoff scales with theta_s (single precision)
        const amrex::Real tol_q =
            (tol * std::abs(m_prof.q())) +
            (10.0_rt * std::numeric_limits<amrex::Real>::epsilon() *
             m_prof.theta_s * m_prof.ustar);
        int nchecked = 0;
        for (int lev = 0; lev < 2; ++lev) {
            amrex::Real tau_along = 0.0_rt;
            amrex::Real tau_cross = 0.0_rt;
            amrex::Real qwall = 0.0_rt;
            amrex::Real ncells = 0.0_rt;
            level_fluxes(lev, tau_along, tau_cross, qwall, ncells);
            if (!(ncells > 0.0_rt)) {
                // Level 0 is covered by level 1 here
                continue;
            }
            ++nchecked;
            EXPECT_NEAR(tau_along, utau2, tol * utau2) << "level " << lev;
            EXPECT_NEAR(tau_cross, 0.0_rt, tol * utau2) << "level " << lev;
            EXPECT_NEAR(qwall, m_prof.q(), tol_q) << "level " << lev;
        }
        ASSERT_EQ(nchecked, 1);
    }

    //! Check that level 0 and level 1 return the same mean stress and heat
    //! flux when level 1 covers part of the wall. The plane averages then
    //! carry the interpolation of the fine-binned profile (2.5 % low at the
    //! level-0 first cell), so both levels differ from the profile's values
    //! by the same factor. The heat flux gets twice the tolerance: the
    //! biased average moves the Obukhov length (by 7 % in the unstable
    //! case), and the two levels feel that through phi_h differently
    void check_levels_agree(const amrex::Real tol)
    {
        ASSERT_NO_FATAL_FAILURE(init_wall_fields());
        amrex::Real tau[2];
        amrex::Real cross[2];
        amrex::Real q[2];
        for (int lev = 0; lev < 2; ++lev) {
            amrex::Real ncells = 0.0_rt;
            level_fluxes(lev, tau[lev], cross[lev], q[lev], ncells);
            ASSERT_GT(ncells, 0.0_rt) << "level " << lev;
        }
        EXPECT_NEAR(tau[1], tau[0], tol * std::abs(tau[0]));
        EXPECT_NEAR(cross[1], cross[0], tol * std::abs(tau[0]));
        EXPECT_NEAR(
            q[1], q[0],
            (2.0_rt * tol * std::abs(m_prof.q())) +
                (10.0_rt * std::numeric_limits<amrex::Real>::epsilon() *
                 m_prof.theta_s * m_prof.ustar));
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
            const amrex::Real wspd = m_prof.speed(z1);
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

    std::string m_refine_box{"0.0 0.0 0.0 60.0 120.0 250.0"};
    bool m_full_coverage{false};
    std::string m_shear_stress_type{"moeng"};
    std::string m_level_means_method;
    bool m_wall_position_set{false};
    amrex::Real m_wall_position{0.0_rt};
    amrex::Real m_log_law_height{0.0_rt};
    //! Specify the surface temperature instead of the surface heat flux
    bool m_surface_temp_mode{false};
    MOProfile m_prof;
    const amrex::Real m_wind_angle{0.5_rt};
    const amrex::Real m_mu{0.01_rt};
};

// With the plane averages read at the exact profile values, every level
// returns u*^2 along the wind and the specified heat flux. With the
// reference-height means on level 1 (ABL.level_means_method = none, the
// behaviour before this change) the moeng stress there is (2 r_m - 1) u*^2
// with r_m = phi_m(3.9 m) / phi_m(11.7 m), about 0.56 u*^2
TEST_F(ABLWallRefinementTest, moeng_wall_model_level_consistent)
{
    m_shear_stress_type = "moeng";
    check_level_fluxes_exact_reference();
}

TEST_F(ABLWallRefinementTest, schumann_wall_model_level_consistent)
{
    m_shear_stress_type = "schumann";
    check_level_fluxes_exact_reference();
}

TEST_F(ABLWallRefinementTest, local_wall_model_level_consistent)
{
    m_shear_stress_type = "local";
    check_level_fluxes_exact_reference();
}

// The constant model uses the plane averages only; this guards that each
// level still returns the friction velocity and the surface heat flux
TEST_F(ABLWallRefinementTest, constant_wall_model_level_consistent)
{
    m_shear_stress_type = "constant";
    check_level_fluxes_exact_reference();
}

// In surface-temperature mode every level keeps the specified surface
// temperature and its mean heat flux is the Monin-Obukhov flux of its own mean
// temperature at its own height; with the reference-height means the level-1
// flux is (r_m + r_h - 1) times that for moeng (0.56 here) and r_h times that
// for schumann and local (0.78), with r_m, r_h the phi_m and phi_h ratios
TEST_F(ABLWallRefinementTest, moeng_surface_temperature_mode_level_consistent)
{
    m_shear_stress_type = "moeng";
    m_surface_temp_mode = true;
    check_level_fluxes_exact_reference();
}

TEST_F(
    ABLWallRefinementTest, schumann_surface_temperature_mode_level_consistent)
{
    m_shear_stress_type = "schumann";
    m_surface_temp_mode = true;
    check_level_fluxes_exact_reference();
}

TEST_F(ABLWallRefinementTest, local_surface_temperature_mode_level_consistent)
{
    m_shear_stress_type = "local";
    m_surface_temp_mode = true;
    check_level_fluxes_exact_reference();
}

// Level 1 covering three quarters of the wall, reference height at the
// level-0 first cell: in neutral conditions the log-law scaling is exact and
// the two levels return the same stress to roundoff (with the
// reference-height means the level-1 stress is 0.69 of the level-0 one for
// moeng and 0.84 for schumann)
TEST_F(ABLWallRefinementTest, levels_agree_at_partial_coverage_neutral)
{
    m_prof.thetastar = 0.0_rt;
    m_shear_stress_type = "moeng";
    check_levels_agree(std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
}

TEST_F(ABLWallRefinementTest, schumann_levels_agree_at_partial_coverage_neutral)
{
    m_prof.thetastar = 0.0_rt;
    m_shear_stress_type = "schumann";
    check_levels_agree(std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
}

// Both levels scaled: with the reference height at the centre of the second
// level-1 cell (11.72 m) and partial coverage, level 0 (first cell 7.8 m) is
// scaled up and level 1 (3.9 m) down along the profile, and in neutral
// conditions the two must still agree to roundoff
TEST_F(ABLWallRefinementTest, levels_agree_with_both_levels_scaled)
{
    m_prof.thetastar = 0.0_rt;
    m_shear_stress_type = "moeng";
    m_log_law_height = 11.71875_rt;
    check_levels_agree(std::numeric_limits<amrex::Real>::epsilon() * 1.0e4_rt);
}

// Unstable conditions with partial coverage: the plane average at the
// reference height is interpolated from the fine-binned profile (2.5 % low),
// so the Obukhov length the wall function finds is 7 % off the profile's and
// the two levels agree to 1 % (stress) and 2 % (heat flux) rather than to
// roundoff; with the reference-height means the level-1 stress is 0.70 of the
// level-0 one
TEST_F(ABLWallRefinementTest, levels_agree_at_partial_coverage_unstable)
{
    m_shear_stress_type = "moeng";
    check_levels_agree(1.0e-2_rt);
}

// ABL.level_means_method = none keeps the reference-height data on every
// level, as before this change; the level-1 stress is then well below the
// level-0 one
TEST_F(ABLWallRefinementTest, none_keeps_reference_height_means)
{
    m_level_means_method = "none";
    m_prof.thetastar = 0.0_rt;
    ASSERT_NO_FATAL_FAILURE(init_wall_fields());
    const auto& abl = sim().physics_manager().get<kynema_sgf::ABL>();
    const auto& wall_func = abl.abl_wall_function();
    EXPECT_EQ(&wall_func.mo(0), &wall_func.mo());
    EXPECT_EQ(&wall_func.mo(1), &wall_func.mo());
    amrex::Real tau[2];
    amrex::Real cross[2];
    amrex::Real q[2];
    for (int lev = 0; lev < 2; ++lev) {
        amrex::Real ncells = 0.0_rt;
        level_fluxes(lev, tau[lev], cross[lev], q[lev], ncells);
        ASSERT_GT(ncells, 0.0_rt) << "level " << lev;
    }
    EXPECT_LT(tau[1], 0.8_rt * tau[0]);
}

TEST_F(ABLWallRefinementTest, invalid_level_means_method_aborts)
{
    m_level_means_method = "own_cells";
    populate_parameters();
    initialize_mesh();
    sim().pde_manager().register_icns();
    sim().create_turbulence_model();
    EXPECT_THROW(sim().init_physics(), amrex::RuntimeError);
}

TEST_F(ABLWallRefinementTest, donelan_wall_model_keeps_reference_height)
{
    // Donelan selects its drag coefficient from the mean wind at the
    // reference height, here 10 m as in the hurricane boundary layer case,
    // with a wind strong enough for the linear range of the coefficient
    m_shear_stress_type = "donelan";
    m_log_law_height = 10.0_rt;
    m_prof.ustar = 1.0_rt;
    check_donelan_reference_height();
}

// A level whose first cell sits at or below the roughness height keeps the
// means at the reference height: the log law, and phi_h, do not hold there
TEST_F(ABLWallRefinementTest, first_cell_below_roughness_keeps_reference)
{
    m_prof.z0 = 5.0_rt; // level-1 first cell at 3.9 m, level-0 at 7.8 m
    m_prof.thetastar = 0.0_rt;
    ASSERT_NO_FATAL_FAILURE(init_wall_fields());
    const auto& abl = sim().physics_manager().get<kynema_sgf::ABL>();
    const auto& wall_func = abl.abl_wall_function();
    ASSERT_NE(&wall_func.mo(1), &wall_func.mo());
    EXPECT_EQ(wall_func.mo(1).zref, wall_func.mo().zref);
    EXPECT_EQ(wall_func.mo(1).vmag_mean, wall_func.mo().vmag_mean);
    EXPECT_GT(wall_func.mo(0).zref, m_prof.z0);
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

// Refinement that does not reach the wall leaves level 0 alone on the wall:
// every level keeps the plane averages at the reference height, as on a
// single-level mesh
TEST_F(ABLWallRefinementTest, refinement_aloft_keeps_reference_height)
{
    m_refine_box = "0.0 0.0 500.0 120.0 120.0 750.0";
    ASSERT_NO_FATAL_FAILURE(init_wall_fields());
    const auto& abl = sim().physics_manager().get<kynema_sgf::ABL>();
    const auto& wall_func = abl.abl_wall_function();
    EXPECT_EQ(&wall_func.mo(0), &wall_func.mo());
    EXPECT_EQ(&wall_func.mo(1), &wall_func.mo());
}

} // namespace kynema_sgf_tests
