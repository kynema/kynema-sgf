#include <AMReX_BC_TYPES.H>
#include "amr-wind/incflo.H"
#include <cmath>

#include "amr-wind/core/Physics.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/physics/BoussinesqBubble.H"
#include "amr-wind/utilities/tagging/RefinementCriteria.H"
#include "amr-wind/utilities/tagging/CartBoxRefinement.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

void incflo::CheckAndSetUpDryRun()
{
    // Check if dry run is requested; exit if not
    {
        amrex::ParmParse pp("incflo");
        pp.query("dry_run", m_dry_run);
        if (!m_dry_run) {
            return;
        }
    }
    // Disable additional computations associated with initialization
    {
        amrex::ParmParse pp("incflo");
        pp.add("initial_iterations", 0);
        pp.add("do_initial_proj", false);
    }
    // Zero time steps, write plotfile and not checkpoint
    {
        amrex::ParmParse pp("time");
        pp.add("max_step", 0);
        pp.add("plot_interval", 1);
        pp.add("checkpoint_inteval", (-1));
    }
    // Give prefix to plotfile
    {
        amrex::ParmParse pp("io");
        std::string current_plt{"plt"};
        pp.query("plot_file", current_plt);
        pp.add("plot_file", (std::string) "dry_run_" + current_plt);
    }
}

/** Parse the input file and populate parameters
 */
void incflo::ReadParameters()
{

    { // Prefix incflo
        amrex::ParmParse pp("incflo");

        pp.query("verbose", m_verbose);

        pp.query("initial_iterations", m_initial_iterations);
        pp.query("do_initial_proj", m_do_initial_proj);

        // Godunov-related flags
        pp.query("use_godunov", m_use_godunov);

        // The default for diffusion_type is 1, i.e. the default m_diff_type is
        // DiffusionType::Crank_Nicolson
        int diffusion_type = 1;
        pp.query("diffusion_type", diffusion_type);
        if (diffusion_type == 0) {
            m_diff_type = DiffusionType::Explicit;
        } else if (diffusion_type == 1) {
            m_diff_type = DiffusionType::Crank_Nicolson;
        } else if (diffusion_type == 2) {
            m_diff_type = DiffusionType::Implicit;
        } else {
            amrex::Abort(
                "We currently require diffusion_type = 0 for explicit, 1 for "
                "Crank-Nicolson or 2 for implicit");
        }

        if (!m_use_godunov && m_time.max_cfl() > 0.5_rt) {
            amrex::Abort(
                "We currently require cfl <= 0.5_rt when using the MOL "
                "advection "
                "scheme");
        }
        if (m_use_godunov && m_time.max_cfl() > 1.0_rt) {
            amrex::Abort(
                "We currently require cfl <= 1.0_rt when using the Godunov "
                "advection scheme");
        }

    } // end prefix incflo

    {
        amrex::ParmParse pp("ICNS");
        pp.query("reconstruct_true_pressure", m_reconstruct_true_pressure);
    }
}

/** Perform initial pressure iterations
 *
 *  Performs a user-defined number of iterations to compute the pressure based
 *  on the initial conditions. This method is only invoked for new simulations
 *  and skipped for restarted simulations from a checkpoint file.
 */
void incflo::InitialIterations()
{
    BL_PROFILE("amr-wind::incflo::InitialIterations()");
    amrex::Print() << "Begin initial pressure iterations. Num. iters = "
                   << m_initial_iterations << '\n';

    compute_dt();

    {
        auto& vel = icns().fields().field;
        vel.copy_state(amr_wind::FieldState::Old, amr_wind::FieldState::New);
        vel.state(amr_wind::FieldState::Old).fillpatch(m_time.current_time());

        if (m_sim.pde_manager().constant_density()) {
            auto& rho = density();
            rho.copy_state(
                amr_wind::FieldState::Old, amr_wind::FieldState::New);
            rho.state(amr_wind::FieldState::Old)
                .fillpatch(m_time.current_time());
        }

        for (auto& eqn : scalar_eqns()) {
            auto& scal = eqn->fields().field;
            scal.copy_state(
                amr_wind::FieldState::Old, amr_wind::FieldState::New);
            scal.state(amr_wind::FieldState::Old)
                .fillpatch(m_time.current_time());
        }
    }

    // Subtract initial pressure if initialized
    // (primarily intended for two-phase overset)
    if (m_reconstruct_true_pressure) {
        auto& press = m_repo.get_field("p");
        amrex::Real global_max_p =
            amr_wind::field_ops::global_max_magnitude(press);
        // If max pressure is different from 0, then infer that pressure has
        // been initialized
        bool init_p = std::abs(global_max_p) > amr_wind::constants::LOOSE_TOL;
        if (init_p) {
            const auto& p0 = m_repo.get_field("reference_pressure");
            for (int lev = 0; lev <= finest_level; lev++) {
                amrex::MultiFab::Subtract(
                    press(lev), p0(lev), 0, 0, 1, p0.num_grow()[0]);
            }
        }
    }

    for (int iter = 0; iter < m_initial_iterations; ++iter) {
        if (m_verbose != 0) {
            amrex::Print() << "In initial_iterations: iter = " << iter << "\n";
        }

        ApplyPredictor(true, 0);

        {
            auto& vel = icns().fields().field;
            // ensure velocity is in stretched mesh space
            if (vel.in_uniform_space() && m_sim.has_mesh_mapping()) {
                vel.to_stretched_space();
            }
            vel.copy_state(
                amr_wind::FieldState::New, amr_wind::FieldState::Old);

            if (m_sim.pde_manager().constant_density()) {
                auto& rho = density();
                rho.copy_state(
                    amr_wind::FieldState::New, amr_wind::FieldState::Old);
            }

            for (auto& eqn : scalar_eqns()) {
                auto& scal = eqn->fields().field;
                scal.copy_state(
                    amr_wind::FieldState::New, amr_wind::FieldState::Old);
            }
        }
    }

    // Add mean pressure back if available
    if (m_reconstruct_true_pressure) {
        auto& press = m_repo.get_field("p");
        const auto& p0 = m_repo.get_field("reference_pressure");
        for (int lev = 0; lev <= finest_level; lev++) {
            amrex::MultiFab::Add(
                press(lev), p0(lev), 0, 0, 1, p0.num_grow()[0]);
        }
    }
    amrex::Print() << "Completed initial pressure iterations" << '\n' << '\n';
}

/** Ensure initial velocity field is divergence-free.
 *
 *  Performs a \ref incflo::ApplyProjection "projection" step to ensure that the
 *  user-provided initial velocity field is divergence free. This method is only
 *  invoked for new simulations. For restarted simulations using a checkpoint
 *  file, this is not necessary.
 */
void incflo::InitialProjection()
{
    BL_PROFILE("amr-wind::incflo::InitialProjection()");

    amrex::Print() << "Begin initial projection" << '\n';
    if (m_verbose != 0) {
        PrintMaxValues("before initial projection");
    }

    amrex::Real dummy_dt = 1.0_rt;
    bool incremental = false;
    ApplyProjection(
        density().vec_const_ptrs(), m_time.current_time(), dummy_dt,
        incremental);

    // We set p and gp back to zero (p0 may still be still non-zero)
    pressure().setVal(0.0_rt);
    grad_p().setVal(0.0_rt);

    if (m_verbose != 0) {
        PrintMaxValues("after initial projection");
    }
    amrex::Print() << "Completed initial projection" << '\n' << '\n';
}
