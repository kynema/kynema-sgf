#include "ks_test_utils/MeshTest.H"
#include "ks_test_utils/iter_tools.H"
#include "ks_test_utils/test_utils.H"
#include "src/equation_systems/vof/vof.H"
#include "src/equation_systems/vof/SplitAdvection.H"
#include "src/equation_systems/vof/split_advection.H"
#include "src/physics/multiphase/MultiPhase.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf_tests {

namespace {
void initialize_volume_fractions(
    const amrex::Box& bx, const amrex::Array4<amrex::Real>& vof_arr)
{
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        vof_arr(i, j, k) = 1.0_rt - (0.1_rt * (i + j + k));
    });
}

void initialize_boundary_volume_fractions(
    const amrex::Box& ghost_box,
    const amrex::Array4<amrex::Real>& values,
    const int poison)
{
    amrex::ParallelFor(ghost_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const int outside =
            (i < 0 || i > 3) + (j < 0 || j > 3) + (k < 0 || k > 3);
        const bool valid = outside <= 1 && i >= -1 && i <= 4 && j >= -1 &&
                           j <= 4 && k >= -1 && k <= 4;
        values(i, j, k) =
            valid ? 0.35_rt + (0.02_rt * i) + (0.04_rt * j) + (0.06_rt * k)
                  : static_cast<amrex::Real>(poison);
    });
}

void initialize_mixed_normal_volume_fractions(
    const amrex::Box& ghost_box, const amrex::Array4<amrex::Real>& values)
{
    amrex::ParallelFor(ghost_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        values(i, j, k) =
            0.4_rt + (0.01_rt * i * i) + (0.02_rt * j) + (0.03_rt * k);
    });
}

void compute_mixed_normal_fluxes(
    const amrex::Box& bx,
    const int dir,
    const int sign,
    const amrex::Array4<amrex::Real const>& fractions,
    const amrex::Array4<amrex::Real>& expected_left,
    const amrex::Array4<amrex::Real>& expected_right)
{
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        amrex::Real mx, my, mz, alpha;
        kynema_sgf::multiphase::fit_plane(
            i, j, k, fractions, mx, my, mz, alpha);
        kynema_sgf::multiphase::eulerian_implicit(
            i, j, k, dir, 0.2_rt, sign, sign, mx, my, mz, alpha, fractions,
            expected_left, expected_right);
    });
}

void reconstruct_neumann_volumes(
    const amrex::Box& cells,
    const amrex::Array4<amrex::Real const>& fractions,
    const amrex::Array4<amrex::Real>& volumes)
{
    amrex::ParallelFor(cells, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        amrex::Real mx, my, mz, alpha;
        kynema_sgf::multiphase::fit_plane(
            i, j, k, i, j, k, fractions, mx, my, mz, alpha);
        volumes(i, j, k) = kynema_sgf::multiphase::cut_volume(
            mx, my, mz, alpha, 0.0_rt, 1.0_rt);
    });
}

void get_accuracy_vofsol(
    kynema_sgf::ScratchField& err_fld,
    const int lev,
    const int dir,
    const int vof_distr,
    const amrex::Real vof_bdyval,
    const kynema_sgf::Field& vof,
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
    kynema_sgf::ScratchField& err_fld,
    const int lev,
    const int dir,
    const int vof_distr,
    const amrex::Real vof_bdyval,
    const kynema_sgf::Field& density,
    const amrex::Real rho1,
    const amrex::Real rho2)
{
    get_accuracy_vofsol(
        err_fld, lev, dir, vof_distr, vof_bdyval, density, rho1, rho2);
}

void get_accuracy_advalpha(
    kynema_sgf::ScratchField& err_fld,
    const int lev,
    const int dir,
    const bool nonzero_flux,
    const amrex::Real rho1,
    const amrex::Real rho2,
    const kynema_sgf::Field& vof,
    const kynema_sgf::Field& advalpha_f)
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
        sim().create_turbulence_model();
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
        auto& multiphase =
            sim().physics_manager().get<kynema_sgf::MultiPhase>();
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
                kynema_sgf::pde::VOF::pde_name() + "-" +
                kynema_sgf::fvm::Godunov::scheme_name());
            seqn.initialize();
            seqn.compute_advection_term(kynema_sgf::FieldState::Old);
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

TEST(VOFBoundaryFlux, inflow_ignores_unfilled_ghosts)
{
    const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(3));
    const amrex::RealBox physical_box(
        {0.0_rt, 0.0_rt, 0.0_rt}, {1.0_rt, 1.0_rt, 1.0_rt});
    const amrex::Array<int, AMREX_SPACEDIM> periodic = {0, 0, 0};
    const amrex::Vector<amrex::Geometry> geom{
        amrex::Geometry(domain, &physical_box, 0, periodic.data())};
    const auto ghost_box = amrex::grow(domain, 2);
    amrex::FArrayBox vof(ghost_box, 1);
    amrex::FArrayBox temporary(amrex::grow(domain, 1), 2);

    for (const auto bc :
         {BC::mass_inflow, BC::mass_inflow_outflow, BC::wave_generation}) {
        amrex::GpuArray<BC, AMREX_SPACEDIM * 2> bcs;
        bcs.fill(bc);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            const auto faces = amrex::surroundingNodes(domain, dir);
            amrex::FArrayBox velocity(ghost_box, 1);
            amrex::FArrayBox flux(faces, 1);
            amrex::FArrayBox advalpha(faces, 1);
            for (const int sign : {-1, 1}) {
                SCOPED_TRACE(
                    ::testing::Message() << "bc=" << static_cast<int>(bc)
                                         << " dir=" << dir << " sign=" << sign);
                velocity.setVal<amrex::RunOn::Device>(
                    static_cast<amrex::Real>(sign));
                amrex::Vector<amrex::Real> reference(faces.numPts());
                amrex::Vector<amrex::Real> result(faces.numPts());
                for (const int poison : {0, 1}) {
                    initialize_boundary_volume_fractions(
                        ghost_box, vof.array(), poison);
                    temporary.setVal<amrex::RunOn::Device>(0.0_rt);
                    kynema_sgf::multiphase::split_compute_fluxes(
                        0, domain, 2 - dir, vof.const_array(),
                        velocity.const_array(), velocity.const_array(),
                        velocity.const_array(), advalpha.array(),
                        advalpha.array(), advalpha.array(), flux.array(),
                        flux.array(), flux.array(), bcs, temporary.dataPtr(),
                        geom, 0.05_rt);
                    auto& output = poison == 0 ? reference : result;
                    amrex::Gpu::copy(
                        amrex::Gpu::deviceToHost, flux.dataPtr(),
                        flux.dataPtr() + faces.numPts(), output.begin());
                }
                amrex::LoopOnCpu(faces, [&](int i, int j, int k) {
                    const amrex::IntVect cell(i, j, k);
                    const int incoming_face = sign > 0 ? 0 : 4;
                    if (cell[dir] == incoming_face) {
                        const auto face = faces.index(cell);
                        EXPECT_NEAR(reference[face], result[face], tol1)
                            << "face=" << face;
                    }
                });
            }
        }
    }
}

TEST(VOFBoundaryFlux, interior_and_periodic_sweeps_use_mixed_normal)
{
    const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(3));
    const amrex::RealBox physical_box(
        {0.0_rt, 0.0_rt, 0.0_rt}, {1.0_rt, 1.0_rt, 1.0_rt});
    const auto ghost_box = amrex::grow(domain, 2);
    const auto work_box = amrex::grow(domain, 1);
    amrex::FArrayBox vof(ghost_box, 1);
    amrex::FArrayBox velocity(ghost_box, 1);
    amrex::FArrayBox actual(work_box, 2);
    amrex::FArrayBox expected(work_box, 2);
    initialize_mixed_normal_volume_fractions(ghost_box, vof.array());
    for (const int is_periodic : {0, 1}) {
        const amrex::Array<int, AMREX_SPACEDIM> periodic = {
            is_periodic, is_periodic, is_periodic};
        const amrex::Geometry geom(domain, &physical_box, 0, periodic.data());
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            for (const int sign : {-1, 1}) {
                SCOPED_TRACE(
                    ::testing::Message() << "periodic=" << is_periodic
                                         << " dir=" << dir << " sign=" << sign);
                velocity.setVal<amrex::RunOn::Device>(
                    static_cast<amrex::Real>(sign));
                actual.setVal<amrex::RunOn::Device>(0.0_rt);
                expected.setVal<amrex::RunOn::Device>(0.0_rt);
                const auto actual_left = actual.array(0);
                const auto actual_right = actual.array(1);
                kynema_sgf::multiphase::sweep_fluxes(
                    dir, domain, 0.2_rt, velocity.const_array(),
                    vof.const_array(), actual_left, actual_right, geom);
                compute_mixed_normal_fluxes(
                    amrex::grow(domain, dir, 1), dir, sign, vof.const_array(),
                    expected.array(0), expected.array(1));
                amrex::Vector<amrex::Real> result(actual.size());
                amrex::Vector<amrex::Real> reference(expected.size());
                amrex::Gpu::copy(
                    amrex::Gpu::deviceToHost, actual.dataPtr(),
                    actual.dataPtr() + actual.size(), result.begin());
                amrex::Gpu::copy(
                    amrex::Gpu::deviceToHost, expected.dataPtr(),
                    expected.dataPtr() + expected.size(), reference.begin());
                const auto check_box = static_cast<bool>(is_periodic)
                                           ? amrex::grow(domain, dir, 1)
                                           : domain;
                amrex::LoopOnCpu(check_box, [&](int i, int j, int k) {
                    const auto cell = work_box.index(amrex::IntVect(i, j, k));
                    for (int component = 0; component < 2; ++component) {
                        const auto index =
                            cell + (component * work_box.numPts());
                        EXPECT_NEAR(result[index], reference[index], tol1)
                            << "cell=" << amrex::IntVect(i, j, k)
                            << " component=" << component;
                    }
                });
            }
        }
    }
}

TEST(VOFBoundaryFlux, neumann_plane_preserves_uniform_fraction)
{
    const amrex::Box cells(amrex::IntVect(-1), amrex::IntVect(1));
    amrex::FArrayBox vof(amrex::grow(cells, 1), 1);
    amrex::FArrayBox reconstructed(cells, 1);
    vof.setVal<amrex::RunOn::Device>(0.35_rt);
    reconstruct_neumann_volumes(
        cells, vof.const_array(), reconstructed.array());
    amrex::Vector<amrex::Real> result(reconstructed.size());
    amrex::Gpu::copy(
        amrex::Gpu::deviceToHost, reconstructed.dataPtr(),
        reconstructed.dataPtr() + reconstructed.size(), result.begin());
    for (const auto volume : result) {
        EXPECT_NEAR(volume, 0.35_rt, tol1);
    }
}

} // namespace kynema_sgf_tests
