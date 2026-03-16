#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/equation_systems/vof/vof.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

namespace {

void init_field3(
    amr_wind::Field& fld,
    const amrex::Real in0,
    const amrex::Real in1,
    const amrex::Real in2)
{
    const int nlevels = fld.repo().num_active_levels();
    const amrex::GpuArray<amrex::Real, 3> ins = {in0, in1, in2};

    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& farrs = fld(lev).arrays();
        amrex::ParallelFor(
            fld(lev), fld.num_grow(), fld.num_comp(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                farrs[nbx](i, j, k, n) = ins[n];
            });
    }
    amrex::Gpu::streamSynchronize();
}

void initialize_volume_fractions(
    const int dir,
    const amrex::Box& bx,
    const amrex::Array4<amrex::Real>& vof_arr)
{

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const int icheck = (dir == 0) ? i : ((dir == 1) ? j : k);
        if (icheck > 0) {
            vof_arr(i, j, k) = 0.0_rt;
        } else {
            vof_arr(i, j, k) = 1.0_rt;
        }
    });
    // Left half is liquid, right half is gas
}

void get_accuracy(
    amrex::Array4<amrex::Real>& err_arr,
    const int dir,
    const amrex::Real rho1,
    const amrex::Real rho2,
    const amrex::Real vel_val,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& problo,
    const amrex::GpuArray<amrex::Real, 3>& varr,
    const amrex::Box& bx,
    const amrex::Array4<const amrex::Real>& um,
    const amrex::Array4<const amrex::Real>& vm,
    const amrex::Array4<const amrex::Real>& wm,
    const amrex::Array4<const amrex::Real>& rf,
    const amrex::Array4<const amrex::Real>& vel,
    const amrex::Array4<const amrex::Real>& dqdt)
{
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const int icheck = (dir == 0) ? i : ((dir == 1) ? j : k);
        // x is face location
        const amrex::Real x = problo[dir] + (icheck * dx[dir]);
        // Check that MAC velocity is as expected, unchanged
        err_arr(i, j, k, 0) = std::abs(um(i, j, k) - varr[0]);
        err_arr(i, j, k, 1) = std::abs(vm(i, j, k) - varr[1]);
        err_arr(i, j, k, 2) = std::abs(wm(i, j, k) - varr[2]);
        // Check that velocity is unchanged after advection
        err_arr(i, j, k, 3) = std::abs(vel(i, j, k, 0) - varr[0]);
        err_arr(i, j, k, 4) = std::abs(vel(i, j, k, 1) - varr[1]);
        err_arr(i, j, k, 5) = std::abs(vel(i, j, k, 2) - varr[2]);

        // Test volume fractions at faces
        if (x == 0.5_rt) {
            // Center face (coming from left cell)
            amrex::Real advvof = 1.0_rt;
            amrex::Real advrho = (rho1 * advvof) + (rho2 * (1.0_rt - advvof));
            err_arr(i, j, k, 6) = std::abs(rf(i, j, k) - advrho);
        } else {
            if (x == 0.0_rt) {
                // Left face (coming from right cell, periodic
                // BC)
                amrex::Real advvof = 0.0_rt;
                amrex::Real advrho =
                    (rho1 * advvof) + (rho2 * (1.0_rt - advvof));
                err_arr(i, j, k, 6) = std::abs(rf(i, j, k) - advrho);
            }
        }

        // Test momentum fluxes by checking convective term
        if (icheck == 0) {
            // Left cell (gas entering, liquid leaving)
            err_arr(i, j, k, 7) = std::abs(
                dqdt(i, j, k, dir) -
                (vel_val * vel_val * (rho2 - rho1) / 0.5_rt));
        } else {
            if (icheck == 1) {
                // Right cell (liquid entering, gas leaving)
                err_arr(i, j, k, 7) = std::abs(
                    dqdt(i, j, k, dir) -
                    (vel_val * vel_val * (rho1 - rho2) / 0.5_rt));
            }
        }
    });
}

} // namespace

class MassMomFluxOpTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{2, 2, 2}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", 2);
            pp.addarr("n_cell", ncell);
        }
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<amrex::Real> problo{{0.0_rt, 0.0_rt, 0.0_rt}};
            amrex::Vector<amrex::Real> probhi{{1.0_rt, 1.0_rt, 1.0_rt}};

            pp.addarr("prob_lo", problo);
            pp.addarr("prob_hi", probhi);
        }
        {
            amrex::ParmParse pp("MultiPhase");
            pp.add("density_fluid1", m_rho1);
            pp.add("density_fluid2", m_rho2);
        }
        {
            amrex::ParmParse pp("incflo");
            amrex::Vector<std::string> physics{"MultiPhase"};
            pp.addarr("physics", physics);
            pp.add("use_godunov", 1);
            pp.add("godunov_type", (std::string) "weno");
        }
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", m_dt);
        }
        {
            amrex::ParmParse pp("transport");
            pp.add("model", (std::string) "TwoPhaseTransport");
            pp.add("viscosity_fluid1", (amrex::Real)0.0_rt);
            pp.add("viscosity_fluid2", (amrex::Real)0.0_rt);
        }
    }

    void testing_coorddir(const int dir)
    {
        constexpr amrex::Real tol =
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e1_rt;

        populate_parameters();
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<int> periodic{{1, 1, 1}};
            pp.addarr("is_periodic", periodic);
        }

        initialize_mesh();

        auto& repo = sim().repo();

        // Set up icns PDE, access to VOF PDE
        auto& pde_mgr = sim().pde_manager();
        auto& mom_eqn = pde_mgr.register_icns();
        mom_eqn.initialize();

        // Initialize physics for the sake of MultiPhase routines
        sim().init_physics();

        // Initialize constant velocity field
        amrex::GpuArray<amrex::Real, 3> varr = {0};
        varr[dir] = m_vel;
        auto& velocity = mom_eqn.fields().field;
        init_field3(velocity, varr[0], varr[1], varr[2]);

        // Initialize volume fraction field
        auto& vof = repo.get_field("vof");
        run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
            auto vof_arr = vof(lev).array(mfi);
            const auto& bx = mfi.validbox();
            initialize_volume_fractions(dir, bx, vof_arr);
        });
        // Populate boundary cells
        vof.fillpatch(0.0_rt);
        // Sync density field with vof
        auto& mphase = sim().physics_manager().get<amr_wind::MultiPhase>();
        mphase.set_density_via_vof();

        // Advance states (new -> old)
        sim().pde_manager().advance_states();

        // Perform pre-advection step to get MAC velocity field (should be
        // uniform)
        mom_eqn.pre_advection_actions(amr_wind::FieldState::Old);

        // Perform VOF solve
        // Get equation handle and perform init
        auto& seqn = pde_mgr(
            amr_wind::pde::VOF::pde_name() + "-" +
            amr_wind::fvm::Godunov::scheme_name());
        seqn.initialize();
        seqn.compute_advection_term(amr_wind::FieldState::Old);
        seqn.post_solve_actions();

        // Zero unused momentum terms: src (pressure)
        auto& grad_p = repo.get_field("gp");
        grad_p.setVal(0.0_rt);

        // Setup mask_cell array to avoid errors in solve
        auto& mask_cell = repo.declare_int_field("mask_cell", 1, 1);
        mask_cell.setVal(1);

        // Perform momentum solve
        mom_eqn.compute_advection_term(amr_wind::FieldState::Old);
        mom_eqn.compute_diffusion_term(amr_wind::FieldState::New);
        mom_eqn.compute_source_term(amr_wind::FieldState::New);
        mom_eqn.compute_predictor_rhs(DiffusionType::Explicit);

        // Get MAC velocities for use in testing
        const auto& umac = repo.get_field("u_mac");
        const auto& vmac = repo.get_field("v_mac");
        const auto& wmac = repo.get_field("w_mac");

        // Get advected alpha
        const std::string sdir = (dir == 0) ? "x" : ((dir == 1) ? "y" : "z");
        const auto& advrho_f = repo.get_field("advalpha_" + sdir);

        // Get convective term
        const auto& conv_term =
            mom_eqn.fields().conv_term.state(amr_wind::FieldState::New);

        // Create scratch field to store error, 1 component for each err type
        auto error_ptr = repo.create_scratch_field(8, 0);
        auto& error_fld = *error_ptr;
        // Initialize at 0
        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            error_fld(lev).setVal(0.0_rt);
        }

        const auto& geom = repo.mesh().Geom();
        run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
            const auto& dx = geom[lev].CellSizeArray();
            const auto& problo = geom[lev].ProbLoArray();
            const auto& bx = mfi.validbox();
            const auto& um = umac(lev).const_array(mfi);
            const auto& vm = vmac(lev).const_array(mfi);
            const auto& wm = wmac(lev).const_array(mfi);
            const auto& rf = advrho_f(lev).const_array(mfi);
            const auto& vel = velocity(lev).const_array(mfi);
            const auto& dqdt = conv_term(lev).const_array(mfi);
            // Populate error scratch field
            auto error = error_fld(lev).array(mfi);
            get_accuracy(
                error, dir, m_rho1, m_rho2, m_vel, dx, problo, varr, bx, um, vm,
                wm, rf, vel, dqdt);
        });

        // Label 0's with variable names so gtest output is decipherable
        constexpr amrex::Real umac_check = 0.0_rt;
        constexpr amrex::Real vmac_check = 0.0_rt;
        constexpr amrex::Real wmac_check = 0.0_rt;
        constexpr amrex::Real u_check = 0.0_rt;
        constexpr amrex::Real v_check = 0.0_rt;
        constexpr amrex::Real w_check = 0.0_rt;
        constexpr amrex::Real advrho_check = 0.0_rt;
        constexpr amrex::Real dqdt_check = 0.0_rt;

        // Check error in each mfab
        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            // Sum error for each component individually and check
            EXPECT_NEAR(error_fld(lev).max(0), umac_check, tol);
            EXPECT_NEAR(error_fld(lev).max(1), vmac_check, tol);
            EXPECT_NEAR(error_fld(lev).max(2), wmac_check, tol);
            EXPECT_NEAR(error_fld(lev).max(3), u_check, tol);
            EXPECT_NEAR(error_fld(lev).max(4), v_check, tol);
            EXPECT_NEAR(error_fld(lev).max(5), w_check, tol);
            EXPECT_NEAR(error_fld(lev).max(6), advrho_check, tol);
            EXPECT_NEAR(error_fld(lev).max(7), dqdt_check, tol);
        }
    }

    const amrex::Real m_rho1 = 1000.0_rt;
    const amrex::Real m_rho2 = 1.0_rt;
    const amrex::Real m_vel = 5.0_rt;
    const amrex::Real m_dt = 0.45_rt * 0.5_rt / m_vel; // first number is CFL
};

TEST_F(MassMomFluxOpTest, fluxfaceX) { testing_coorddir(0); }
TEST_F(MassMomFluxOpTest, fluxfaceY) { testing_coorddir(1); }
TEST_F(MassMomFluxOpTest, fluxfaceZ) { testing_coorddir(2); }

} // namespace amr_wind_tests
