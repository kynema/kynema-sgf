#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/equation_systems/vof/vof.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {

namespace {
void initialize_volume_fractions(
    const amrex::Box& bx, const amrex::Array4<amrex::Real>& vof_arr)
{
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        vof_arr(i, j, k) = 1.0_rt - (0.1_rt * (i + j + k));
    });
}

void get_accuracy_vofsol(
    amr_wind::ScratchField& err_fld,
    const int lev,
    const int dir,
    const int vof_distr,
    const amrex::Real vof_bdyval,
    const amr_wind::Field& vof,
    const amrex::Real wt1 = 1.0_rt,
    const amrex::Real wt2 = 0.0_rt)
{
    /* -- Check VOF boundary values from fillpatch -- */
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(vof(lev), amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
        auto err_arr = err_fld(lev).array(mfi);
        const auto& vof_arr = vof(lev).array(mfi);
        // Check lo and hi
        amrex::ParallelFor(2, [=] AMREX_GPU_DEVICE(int i0) {
            int off = i0 * 3;
            int i = 0;
            int j = 0;
            int k = 0;
            // Small mesh, loop in serial for check
            for (int i1 = 0; i1 < 2; ++i1) {
                for (int i2 = 0; i2 < 2; ++i2) {
                    // Organize indices
                    if (dir == 0) {
                        i = -1 + off;
                        j = i1;
                        k = i2;
                    } else if (dir == 1) {
                        i = i1;
                        j = -1 + off;
                        k = i2;
                    } else if (dir == 2) {
                        i = i1;
                        j = i2;
                        k = -1 + off;
                    }
                    // Calculate reference value
                    amrex::Real ref_val = vof_bdyval; // case 0
                    if (vof_distr == 1) {
                        // hoextrap
                        ref_val = 1.0_rt - (0.1_rt * (i + j + k));
                    } else if (vof_distr == 2) {
                        // foextrap
                        const int ii = amrex::max(0, amrex::min(1, i));
                        const int jj = amrex::max(0, amrex::min(1, j));
                        const int kk = amrex::max(0, amrex::min(1, k));
                        ref_val = 1.0_rt - (0.1_rt * (ii + jj + kk));
                    }
                    // Check against reference value
                    const amrex::Real wt_ref =
                        (wt1 * ref_val) + (wt2 * (1.0_rt - ref_val));
                    err_arr(i, j, k, 0) = std::abs(vof_arr(i, j, k) - wt_ref);
                }
            }
        });
    }
}

void get_accuracy_density(
    amr_wind::ScratchField& err_fld,
    const int lev,
    const int dir,
    const int vof_distr,
    const amrex::Real vof_bdyval,
    const amr_wind::Field& density,
    const amrex::Real rho1,
    const amrex::Real rho2)
{
    get_accuracy_vofsol(
        err_fld, lev, dir, vof_distr, vof_bdyval, density, rho1, rho2);
}

void get_accuracy_advalpha(
    amr_wind::ScratchField& err_fld,
    const int lev,
    const int dir,
    const bool nonzero_flux,
    const amrex::Real rho1,
    const amrex::Real rho2,
    const amr_wind::Field& vof,
    const amr_wind::Field& advalpha_f)
{
    /* -- Check VOF boundary fluxes -- */
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(vof(lev), amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

        auto err_arr = err_fld(lev).array(mfi);
        const auto& af = advalpha_f(lev).array(mfi);
        // Check lo and hi
        amrex::ParallelFor(2, [=] AMREX_GPU_DEVICE(int i0) {
            int off = i0 * 2;
            int i = 0;
            int j = 0;
            int k = 0;
            // Small mesh, loop in serial for check
            for (int i1 = 0; i1 < 2; ++i1) {
                for (int i2 = 0; i2 < 2; ++i2) {
                    // Organize indices
                    if (dir == 0) {
                        i = 0 + off;
                        j = i1;
                        k = i2;
                    } else if (dir == 1) {
                        i = i1;
                        j = 0 + off;
                        k = i2;
                    } else if (dir == 2) {
                        i = i1;
                        j = i2;
                        k = 0 + off;
                    }
                    // Check whether flux is nonzero
                    constexpr amrex::Real advvof = 0.0_rt;
                    const amrex::Real advrho =
                        (rho1 * advvof) + (rho2 * (1.0_rt - advvof));
                    if (nonzero_flux) {
                        err_arr(i, j, k, 0) =
                            af(i, j, k) > 0.0_rt ? 0.0_rt : 1.0_rt;
                    } else {
                        err_arr(i, j, k, 0) = std::abs(af(i, j, k) - advrho);
                    }
                }
            }
        });
    }
}
} // namespace

class VOFBCTest : public MeshTest
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
    }

    void
    testing_bc_coorddir(const int option, const int dir, const amrex::Real tol)
    {

        // Get string version of direction
        const std::string sdir = (dir == 0) ? "x" : ((dir == 1) ? "y" : "z");
        const std::string odir0 = (dir == 0) ? "y" : "x";
        const std::string odir1 = (dir == 2) ? "y" : "z";

        // Set up "correct" answers according to each option
        int vof_distr = 0;
        bool nonzero_flux = true;
        std::string bc_string = "sixteen chars --";
        if (option == 1) {
            // Mass inflow
            vof_distr = 0;       // prescribed vof value in bdy
            nonzero_flux = true; // flux can occur
            bc_string = "mass_inflow";
        } else if (option == 2) {
            // Slip wall
            vof_distr = 2;        // first-order extrapolation
            nonzero_flux = false; // flux cannot occur
            bc_string = "slip_wall";
        } else if (option == 3) {
            // No-slip wall
            vof_distr = 2;        // zero-gradient at bdy
            nonzero_flux = false; // flux cannot occur
            bc_string = "no_slip_wall";
        } else if (option == 4) {
            // Pressure outflow
            vof_distr = 2;       // zero-gradient at bdy
            nonzero_flux = true; // flux can occur
            bc_string = "pressure_outflow";
        }

        populate_parameters();
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<int> periodic{{0, 0, 0}};
            pp.addarr("is_periodic", periodic);
        }
        {
            amrex::ParmParse pp(sdir + "lo");
            pp.add("type", bc_string);
            if (option == 1) {
                // Specify vof
                pp.add("vof", m_vof_bdyval);
                pp.add("density", m_rho1);
                // Specify other quantities to avoid errors, not actually used
                pp.add("levelset", 0.0_rt);
                pp.addarr(
                    "velocity",
                    amrex::Vector<amrex::Real>{0.0_rt, 0.0_rt, 0.0_rt});
            }
        }
        {
            amrex::ParmParse pp(sdir + "hi");
            pp.add("type", bc_string);
            if (option == 1) {
                pp.add("vof", m_vof_bdyval);
                pp.add("density", m_rho1);
                pp.add("levelset", 0.0_rt);
                pp.addarr(
                    "velocity",
                    amrex::Vector<amrex::Real>{0.0_rt, 0.0_rt, 0.0_rt});
            }
        }
        {
            amrex::ParmParse pp(odir0 + "lo");
            pp.add("type", (std::string) "slip_wall");
        }
        {
            amrex::ParmParse pp(odir0 + "hi");
            pp.add("type", (std::string) "slip_wall");
        }
        {
            amrex::ParmParse pp(odir1 + "lo");
            pp.add("type", (std::string) "slip_wall");
        }
        {
            amrex::ParmParse pp(odir1 + "hi");
            pp.add("type", (std::string) "slip_wall");
        }

        initialize_mesh();

        auto& repo = sim().repo();

        // Set up icns PDE to set up mac velocity, access to VOF PDE
        auto& pde_mgr = sim().pde_manager();
        auto& mom_eqn = pde_mgr.register_icns();
        mom_eqn.initialize();

        // Initialize physics for the sake of MultiPhase routines
        sim().init_physics();

        // Initialize volume fraction field
        auto& vof = repo.get_field("vof");
        run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
            auto vof_arr = vof(lev).array(mfi);
            const auto& bx = mfi.validbox();
            initialize_volume_fractions(bx, vof_arr);
        });

        // Populate boundary cells
        vof.fillpatch(0.0_rt);

        // Base level
        int lev = 0;

        // Create scratch field to store error (with ghost cells for BCs)
        auto error_ptr = repo.create_scratch_field(1, 1);
        auto& error_fld = *error_ptr;
        // Initialize at 0
        error_fld(lev).setVal(0.0_rt);

        // Compute error
        get_accuracy_vofsol(error_fld, lev, dir, vof_distr, m_vof_bdyval, vof);

        // Check error from first part
        constexpr amrex::Real refval_check = 0.0_rt;
        const amrex::Real vof_distr_err = error_fld(lev).max(0, 1);
        EXPECT_NEAR(vof_distr_err, refval_check, tol);

        // Calculate density and fillpatch
        auto& density = repo.get_field("density");
        auto& multiphase = sim().physics_manager().get<amr_wind::MultiPhase>();
        multiphase.set_density_via_vof();
        density.fillpatch(0.0_rt);

        // Compute error of density field and check
        error_fld(lev).setVal(0.0_rt);
        get_accuracy_density(
            error_fld, lev, dir, vof_distr, m_vof_bdyval, density, m_rho1,
            m_rho2);
        const amrex::Real density_distr_err = error_fld(lev).max(0, 1);
        EXPECT_NEAR(density_distr_err, refval_check, tol);

        // Test positive and negative velocity
        for (int sign = -1; sign < 2; sign += 2) {
            // Reset error for each sign value
            error_fld(lev).setVal(0.0_rt);

            // Get mac velocity fields and set values based on dir
            auto& umac = repo.get_field("u_mac");
            auto& vmac = repo.get_field("v_mac");
            auto& wmac = repo.get_field("w_mac");
            umac.setVal(dir == 0 ? m_vel * sign : 0.0_rt);
            vmac.setVal(dir == 1 ? m_vel * sign : 0.0_rt);
            wmac.setVal(dir == 2 ? m_vel * sign : 0.0_rt);

            // Perform VOF solve
            // Advance states (new -> old)
            sim().pde_manager().advance_states();
            // Get equation handle and perform init
            auto& seqn = pde_mgr(
                amr_wind::pde::VOF::pde_name() + "-" +
                amr_wind::fvm::Godunov::scheme_name());
            seqn.initialize();
            seqn.compute_advection_term(amr_wind::FieldState::Old);
            seqn.post_solve_actions();

            // Get advected alpha
            const auto& advalpha_f = repo.get_field("advalpha_" + sdir);

            // Compute error
            get_accuracy_advalpha(
                error_fld, lev, dir, nonzero_flux, m_rho1, m_rho2, vof,
                advalpha_f);

            // Check error from second part for this value of sign
            constexpr amrex::Real flux_check = 0.0_rt;
            EXPECT_NEAR(error_fld(lev).max(0, 1), flux_check, tol);
        }
    }
    const amrex::Real m_rho1 = 1000.0_rt;
    const amrex::Real m_rho2 = 1.0_rt;
    const amrex::Real m_vel = 5.0_rt;
    const amrex::Real m_vof_bdyval = 1.0_rt;
    const amrex::Real m_dt = 0.45_rt * 0.5_rt / m_vel; // first number is CFL
};

constexpr amrex::Real tol1 =
    std::numeric_limits<amrex::Real>::epsilon() * 1.0e1_rt;

TEST_F(VOFBCTest, dirichletX) { testing_bc_coorddir(1, 0, tol1); }
TEST_F(VOFBCTest, slipwallY) { testing_bc_coorddir(2, 1, tol1); }
TEST_F(VOFBCTest, noslipwallZ) { testing_bc_coorddir(3, 2, tol1); }
TEST_F(VOFBCTest, pressureX) { testing_bc_coorddir(4, 0, tol1); }

} // namespace amr_wind_tests
