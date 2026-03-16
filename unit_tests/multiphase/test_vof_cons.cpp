#include "aw_test_utils/MeshTest.H"
#include "aw_test_utils/iter_tools.H"
#include "aw_test_utils/test_utils.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/equation_systems/vof/vof.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/utilities/tagging/CartBoxRefinement.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind_tests {
namespace {

void initialize_volume_fractions(
    const int dir, const int nx, amr_wind::Field& vof)
{

    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto vof_arr = vof(lev).array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (dir < 0) {
                // Bottom left half is liquid, top right half is gas
                if (i + j + k == 3) {
                    vof_arr(i, j, k) = 0.5_rt;
                } else {
                    if (i + j + k < 3) {
                        vof_arr(i, j, k) = 1.0_rt;
                    } else {
                        vof_arr(i, j, k) = 0.0_rt;
                    }
                }
            } else {
                const int icheck = (dir == 0) ? i : ((dir == 1) ? j : k);
                if ((2 * icheck) + 1 == nx) {
                    vof_arr(i, j, k) = 0.5_rt;
                } else {
                    if ((2 * icheck) + 1 < nx) {
                        vof_arr(i, j, k) = 1.0_rt;
                    } else {
                        vof_arr(i, j, k) = 0.0_rt;
                    }
                }
            }
        });
    });
    // Populate boundary cells
    vof.fillpatch(0.0_rt);
}

void initialize_adv_velocities(
    amr_wind::Field& vof,
    amr_wind::Field& umac,
    amr_wind::Field& vmac,
    amr_wind::Field& wmac,
    amrex::GpuArray<amrex::Real, 3> varr)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto um = umac(lev).array(mfi);
        auto vm = vmac(lev).array(mfi);
        auto wm = wmac(lev).array(mfi);
        const auto& gbx = mfi.growntilebox(1);
        amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            um(i, j, k) = varr[0];
            vm(i, j, k) = varr[1];
            wm(i, j, k) = varr[2];
        });
    });
}

void get_accuracy(
    amr_wind::ScratchField& err_fld, int dir, int nx, amr_wind::Field& vof)
{
    run_algorithm(vof, [&](const int lev, const amrex::MFIter& mfi) {
        auto err_arr = err_fld(lev).array(mfi);
        const auto& vof_arr = vof(lev).const_array(mfi);
        const auto& bx = mfi.validbox();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const int icheck = (dir == 0) ? i : ((dir == 1) ? j : k);
            // Check if current solution matches initial solution
            if ((2 * icheck) + 1 == nx) {
                err_arr(i, j, k) = std::abs(vof_arr(i, j, k) - 0.5_rt);
            } else {
                if ((2 * icheck) + 1 < nx) {
                    err_arr(i, j, k) = std::abs(vof_arr(i, j, k) - 1.0_rt);
                } else {
                    err_arr(i, j, k) = std::abs(vof_arr(i, j, k) - 0.0_rt);
                }
            }
        });
    });
}
} // namespace

class VOFConsTest : public MeshTest
{
protected:
    void populate_parameters() override
    {
        MeshTest::populate_parameters();

        {
            amrex::ParmParse pp("amr");
            amrex::Vector<int> ncell{{m_nx, m_nx, m_nx}};
            pp.add("max_level", 0);
            pp.add("max_grid_size", m_nx);
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
        }
        {
            amrex::ParmParse pp("time");
            pp.add("fixed_dt", m_dt);
        }
        {
            amrex::ParmParse pp("VOF");
            pp.add("remove_debris", 0);
        }
    }

    void testing_coorddir(const int dir, amrex::Real CFL)
    {
        constexpr amrex::Real tol =
            std::numeric_limits<amrex::Real>::epsilon() * 1.0e1_rt;

        // Flow-through time
        const amrex::Real ft_time = 1.0_rt / m_vel;

        // Set timestep according to input
        m_dt = ft_time / ((amrex::Real)m_nx) * CFL;
        // Round to nearest integer timesteps
        int niter = (int)round(ft_time / m_dt);
        // Modify dt to fit niter
        m_dt = ft_time / ((amrex::Real)niter);

        populate_parameters();
        {
            amrex::ParmParse pp("geometry");
            amrex::Vector<int> periodic{{1, 1, 1}};
            pp.addarr("is_periodic", periodic);
        }
        // dir = -2 corresponds to multi-level
        if (dir == -2) {
            {
                amrex::ParmParse pp("amr");
                amrex::Vector<int> ncell{{m_nx + 1, m_nx + 1, m_nx + 1}};
                pp.add("max_level", 1);
                pp.add("max_grid_size", m_nx + 1);
                pp.add("blocking_factor", 2);
                pp.addarr("n_cell", ncell);
            }
            // Create the "input file"
            std::stringstream ss;
            ss << "1 // Number of levels" << '\n';
            ss << "1 // Number of boxes at this level" << '\n';
            ss << "0.8 0.5 0.5 0.9 0.5 0.5" << '\n';

            create_mesh_instance<RefineMesh>();
            std::unique_ptr<amr_wind::CartBoxRefinement> box_refine(
                new amr_wind::CartBoxRefinement(sim()));
            box_refine->read_inputs(mesh(), ss);

            if (mesh<RefineMesh>() != nullptr) {
                mesh<RefineMesh>()->refine_criteria_vec().push_back(
                    std::move(box_refine));
            }
        }

        initialize_mesh();

        auto& repo = sim().repo();

        // PDE manager, for access to VOF PDE later
        auto& pde_mgr = sim().pde_manager();
        // Setup of icns provides MAC velocities
        pde_mgr.register_icns();

        // Initialize physics for the sake of MultiPhase routines
        sim().init_physics();

        // Initialize volume fraction field
        auto& vof = repo.get_field("vof");
        initialize_volume_fractions(dir, m_nx, vof);
        // Get multiphase object
        auto& mphase = sim().physics_manager().get<amr_wind::MultiPhase>();

        // Initialize constant velocity field in single direction
        amrex::GpuArray<amrex::Real, 3> varr = {0};
        if (dir < 0) {
            varr[0] = m_vel;
            varr[1] = m_vel;
            varr[2] = m_vel;
        } else {
            varr[dir] = m_vel;
        }

        // Set advection velocities accordingly
        auto& umac = repo.get_field("u_mac");
        auto& vmac = repo.get_field("v_mac");
        auto& wmac = repo.get_field("w_mac");
        initialize_adv_velocities(vof, umac, vmac, wmac, varr);

        // Get initial VOF sum
        amrex::Real sum_vof0 = mphase.volume_fraction_sum();
        // Get equation handle and perform init
        auto& seqn = pde_mgr(
            amr_wind::pde::VOF::pde_name() + "-" +
            amr_wind::fvm::Godunov::scheme_name());
        seqn.initialize();

        for (int n = 0; n < niter; ++n) {
            // Copy new to old to prep for advection
            for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
                amrex::MultiFab::Copy(
                    vof.state(amr_wind::FieldState::Old)(lev), vof(lev), 0, 0,
                    vof.num_comp(), vof.num_grow());
            }
            // Perform VOF solve
            seqn.compute_advection_term(amr_wind::FieldState::Old);
            seqn.post_solve_actions();
            // Check conservation
            EXPECT_NEAR(mphase.volume_fraction_sum(), sum_vof0, tol);
        }

        if (dir >= 0) {
            // Create scratch field to store error
            auto error_ptr = repo.create_scratch_field(1, 0);
            auto& error_fld = *error_ptr;
            // Initialize at 0
            for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
                error_fld(lev).setVal(0.0_rt);
            }

            get_accuracy(error_fld, dir, m_nx, vof);

            // Check error in each mfab
            constexpr amrex::Real vofsol_check = 0.0_rt;
            for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
                // Sum error and check
                EXPECT_NEAR(error_fld(lev).max(0), vofsol_check, tol);
            }
        }
    }
    const amrex::Real m_rho1 = 1000.0_rt;
    const amrex::Real m_rho2 = 1.0_rt;
    const amrex::Real m_vel = 5.0_rt;
    const int m_nx = 3;
    amrex::Real m_dt = 0.0_rt; // will be set according to CFL
};

TEST_F(VOFConsTest, X) { testing_coorddir(0, 0.45_rt); }
TEST_F(VOFConsTest, Y) { testing_coorddir(1, 0.45_rt); }
TEST_F(VOFConsTest, Z) { testing_coorddir(2, 0.45_rt); }
// Need multi-directional velocity and vof field to test communication of vof
// during directionally-split advection
TEST_F(VOFConsTest, CFL045) { testing_coorddir(-1, 0.45_rt); }
TEST_F(VOFConsTest, CFL01) { testing_coorddir(-1, 0.1_rt); }
// Test transport across multiple mesh levels - just check conservation
TEST_F(VOFConsTest, 2level) { testing_coorddir(-2, 0.5_rt * 0.45_rt); }

} // namespace amr_wind_tests
