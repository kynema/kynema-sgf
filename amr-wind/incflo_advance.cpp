#include <cmath>

#include "amr-wind/incflo.H"
#include "amr-wind/core/Physics.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/equation_systems/PDEBase.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include "amr-wind/utilities/console_io.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

void incflo::pre_advance_stage1()
{
    BL_PROFILE("amr-wind::incflo::pre_advance_stage1");
    advance_time();
}

void incflo::pre_advance_stage2()
{
    BL_PROFILE("amr-wind::incflo::pre_advance_stage2");
    for (auto& pp : m_sim.physics()) {
        pp->pre_advance_work();
    }

    m_sim.helics().pre_advance_work();
}

void incflo::prepare_time_step()
{
    m_sim.pde_manager().advance_states();
    m_sim.pde_manager().prepare_boundaries();
    for (auto& pp : m_sim.physics()) {
        pp->pre_predictor_work();
    }
}

/** Advance simulation state by one timestep
 *
 *  Performs the following actions at a given timestep
 *  - Compute \f$\Delta t\f$
 *  - Advance all computational fields to new timestate in preparation for time
 * integration
 *  - Call pre-advance work for all registered physics modules
 *  - For Godunov scheme, advance to new time state
 *  - For MOL scheme, call predictor corrector steps
 *  - Perform any post-advance work
 *
 *  Much of the heavy-lifting is done by incflo::ApplyPredictor and
 *  incflo::ApplyCorrector. Please refer to the documentation of those methods
 *  for detailed information on the various equations being solved.
 *
 * \callgraph
 */
void incflo::advance(const int fixed_point_iteration)
{
    BL_PROFILE("amr-wind::incflo::advance");
    if (fixed_point_iteration == 0) {
        prepare_time_step();
    }

    ApplyPredictor(false, fixed_point_iteration);

    if (!m_use_godunov) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            fixed_point_iteration == 0,
            "Fixed point iterations are not supported for MOL");

        ApplyCorrector();
    }
}

// Apply predictor step
//
//  For Godunov, this completes the timestep. For MOL, this is the first part of
//  the predictor/corrector within a timestep.
//
//  <ol>
//  <li> Use u = vel_old to compute
//
//     \code{.cpp}
//     conv_u  = - u grad u
//     conv_r  = - div( u rho  )
//     conv_t  = - div( u trac )
//     eta_old     = visosity at m_time.current_time()
//     if (m_diff_type == DiffusionType::Explicit)
//        divtau _old = div( eta ( (grad u) + (grad u)^T ) ) / rho^n
//        rhs = u + dt * ( conv + divtau_old )
//     else
//        divtau_old  = 0.0_rt
//        rhs = u + dt * conv
//
//     eta     = eta at new_time
//     \endcode
//
//  <li> Add explicit forcing term i.e. gravity + lagged pressure gradient
//
//     \code{.cpp}
//     rhs += dt * ( g - grad(p + p0) / rho^nph )
//     \endcode
//
//  Note that in order to add the pressure gradient terms divided by rho,
//  we convert the velocity to momentum before adding and then convert them
//  back.
//
//  <li> A. If (m_diff_type == DiffusionType::Implicit)
//        solve implicit diffusion equation for u*
//
//  \code{.cpp}
//  ( 1 - dt / rho^nph * div ( eta grad ) ) u* = u^n + dt * conv_u
//                                               + dt * ( g - grad(p + p0) /
//                                               rho^nph )
//  \endcode
//
//  B. If (m_diff_type == DiffusionType::Crank-Nicolson)
//     solve semi-implicit diffusion equation for u*
//
//     \code{.cpp}
//     ( 1 - (dt/2) / rho^nph * div ( eta_old grad ) ) u* = u^n +
//            dt * conv_u + (dt/2) / rho * div (eta_old grad) u^n
//          + dt * ( g - grad(p + p0) / rho^nph )
//     \endcode
//
//  <li> Apply projection (see incflo::ApplyProjection)
//
//     Add pressure gradient term back to u*:
//
//      \code{.cpp}
//      u** = u* + dt * grad p / rho^nph
//      \endcode
//
//     Solve Poisson equation for phi:
//
//     \code{.cpp}
//     div( grad(phi) / rho^nph ) = div( u** )
//     \endcode
//
//     Update pressure:
//
//     p = phi / dt
//
//     Update velocity, now divergence free
//
//     vel = u** - dt * grad p / rho^nph
//  </ol>
//
// It is assumed that the ghost cels of the old data have been filled and
// the old and new data are the same in valid region.
//

/** Apply predictor step
 *
 *  For Godunov, this completes the timestep. For MOL, this is the first part of
 *  the predictor/corrector within a timestep.
 *
 *  <ol>
 *  <li> Solve transport equation for momentum and scalars
 *
 *  \f{align}
 *  \left[1 - \kappa \frac{\Delta t}{\rho^{n+1/2}} \nabla \cdot \left( \mu
 *  \nabla \right)\right] u^{*} &= u^n - \Delta t (u \cdot \nabla) u + (1 -
 * \kappa) \frac{\Delta t}{\rho^n} \nabla \cdot \left( \mu^{n} \nabla\right)
 * u^{n} + \frac{\Delta t}{\rho^{n+1/2}} \left( S_u - \nabla(p + p_0)\right) \\
 *  \f}
 *
 *  where
 *  \f{align}
 *  \kappa = \begin{cases}
 *  0 & \text{Explicit} \\
 *  0.5 & \text{Crank-Nicolson} \\
 *  1 & \text{Implicit}
 *  \end{cases}
 *  \f}
 *
 *  <li> \ref incflo::ApplyProjection "Apply projection"
 *  </ol>
 */
void incflo::ApplyPredictor(
    const bool incremental_projection, const int fixed_point_iteration)
{
    BL_PROFILE("amr-wind::incflo::ApplyPredictor");
    // We use the new time value for things computed on the "*" state
    amrex::Real new_time = m_time.new_time();

    if (m_verbose > 2) {
        PrintMaxValues("before predictor step");
    }

    if (m_use_godunov) {
        amr_wind::io::print_mlmg_header("Godunov:");
    } else {
        amr_wind::io::print_mlmg_header("Predictor:");
    }

    auto& icns_fields = icns().fields();
    auto& velocity_new = icns_fields.field;
    auto& velocity_old = velocity_new.state(amr_wind::FieldState::Old);

    auto& density_new = density();
    const auto& density_old = density_new.state(amr_wind::FieldState::Old);
    auto& density_nph = density_new.state(amr_wind::FieldState::NPH);

    if (fixed_point_iteration > 0) {
        icns().fields().field.fillpatch(m_time.current_time());
        // Get n + 1/2 velocity
        amr_wind::field_ops::lincomb(
            icns().fields().field.state(amr_wind::FieldState::NPH), 0.5_rt,
            icns().fields().field.state(amr_wind::FieldState::Old), 0, 0.5_rt,
            icns().fields().field, 0, 0, icns().fields().field.num_comp(),
            icns().fields().field.num_grow());
    }

    // *************************************************************************************
    // Compute viscosity / diffusive coefficients
    // *************************************************************************************
    // TODO: This sub-section has not been adjusted for mesh mapping - adjust in
    // corrector too
    m_sim.turbulence_model().update_turbulent_viscosity(
        amr_wind::FieldState::Old, m_diff_type);
    icns().compute_mueff(amr_wind::FieldState::Old);
    for (auto& eqns : scalar_eqns()) {
        eqns->compute_mueff(amr_wind::FieldState::Old);
    }

    // *************************************************************************************
    // Define the forcing terms to use in the Godunov prediction
    // *************************************************************************************
    // TODO: Godunov has not been adjusted for mesh mapping - adjust in
    // corrector too
    if (m_use_godunov) {
        icns().compute_source_term(amr_wind::FieldState::Old);
        for (auto& seqn : scalar_eqns()) {
            seqn->compute_source_term(amr_wind::FieldState::Old);
        }
    }

    // TODO: This sub-section has not been adjusted for mesh mapping - adjust in
    // corrector too
    if (need_divtau()) {
        // *************************************************************************************
        // Compute explicit viscous term using old density (1/rho)
        // *************************************************************************************
        // Reuse existing buffer to avoid creating new multifabs
        amr_wind::field_ops::copy(
            velocity_new, velocity_old, 0, 0, velocity_new.num_comp(), 1);
        icns().compute_diffusion_term(amr_wind::FieldState::Old);
        if (m_use_godunov) {
            auto& velocity_forces = icns_fields.src_term;
            // only the old states are used in predictor
            const auto& divtau = icns_fields.diff_term;

            amr_wind::field_ops::add(
                velocity_forces, divtau, 0, 0, AMREX_SPACEDIM, 0);
        }
        // *************************************************************************************
        // Compute explicit diffusive terms
        // *************************************************************************************
        for (auto& eqn : scalar_eqns()) {
            auto& field = eqn->fields().field;
            // Reuse existing buffer to avoid creating new multifabs
            amr_wind::field_ops::copy(
                field, field.state(amr_wind::FieldState::Old), 0, 0,
                field.num_comp(), 1);

            eqn->compute_diffusion_term(amr_wind::FieldState::Old);

            if (m_use_godunov) {
                amr_wind::field_ops::add(
                    eqn->fields().src_term, eqn->fields().diff_term, 0, 0,
                    field.num_comp(), 0);
            }
        }
    }

    if (m_use_godunov) {
        // ICNS source term (and viscous term) are in terms of momentum;
        // convert to velocity for MAC velocities by dividing by density
        amr_wind::field_ops::divide(
            icns().fields().src_term, density_old, 0, 0, 1, AMREX_SPACEDIM, 0);

        const int nghost_force = 1;
        amrex::IntVect ng(nghost_force);
        icns().fields().src_term.fillpatch(m_time.current_time(), ng);

        for (auto& eqn : scalar_eqns()) {
            eqn->fields().src_term.fillpatch(m_time.current_time(), ng);
        }
    }

    // Extrapolate and apply MAC projection for advection velocities
    const auto fstate_preadv = (fixed_point_iteration == 0)
                                   ? amr_wind::FieldState::Old
                                   : amr_wind::FieldState::NPH;
    icns().pre_advection_actions(fstate_preadv);

    // *************************************************************************************
    // Update density first
    // *************************************************************************************
    if (m_sim.pde_manager().constant_density()) {
        amr_wind::field_ops::copy(density_nph, density_old, 0, 0, 1, 1);
    }

    // TODO: This sub-section has not been adjusted for mesh mapping - adjust in
    // corrector too.
    // Perform scalar update one at a time. This is to allow an
    // updated density at `n+1/2` to be computed before other scalars use it
    // when computing their source terms.
    for (auto& eqn : scalar_eqns()) {
        // Compute explicit advection
        eqn->compute_advection_term(amr_wind::FieldState::Old);

        // Compute (recompute for Godunov) the scalar forcing terms
        eqn->compute_source_term(amr_wind::FieldState::NPH);

        // Update the scalar (if explicit), or the RHS for implicit/CN
        eqn->compute_predictor_rhs(m_diff_type);

        auto& field = eqn->fields().field;
        if (m_diff_type != DiffusionType::Explicit) {
            amrex::Real dt_diff = (m_diff_type == DiffusionType::Implicit)
                                      ? m_time.delta_t()
                                      : 0.5_rt * m_time.delta_t();

            // Solve diffusion eqn. and update of the scalar field
            eqn->solve(dt_diff);

            // Post-processing actions after a PDE solve
        } else if (m_diff_type == DiffusionType::Explicit && m_use_godunov) {
            // explicit RK2
            auto diff_old =
                m_repo.create_scratch_field(1, 0, amr_wind::FieldLoc::CELL);
            auto& diff_new =
                eqn->fields().diff_term.state(amr_wind::FieldState::New);
            amr_wind::field_ops::copy(*diff_old, diff_new, 0, 0, 1, 0);
            eqn->compute_diffusion_term(amr_wind::FieldState::New);
            amr_wind::field_ops::saxpy(
                diff_new, -1.0_rt, *diff_old, 0, 0, 1, 0);
            eqn->improve_explicit_diffusion(m_time.delta_t());
        }
        eqn->post_solve_actions();

        // Update scalar at n+1/2
        amr_wind::field_ops::lincomb(
            field.state(amr_wind::FieldState::NPH), 0.5_rt,
            field.state(amr_wind::FieldState::Old), 0, 0.5_rt, field, 0, 0,
            field.num_comp(), 1);
    }

    // With scalars computed, compute advection of momentum
    const auto fstate = (fixed_point_iteration == 0)
                            ? amr_wind::FieldState::Old
                            : amr_wind::FieldState::NPH;
    icns().compute_advection_term(fstate);

    // *************************************************************************************
    // Define (or if use_godunov, re-define) the forcing terms and viscous
    // terms independently for the right hand side, without 1/rho term
    // *************************************************************************************
    icns().compute_source_term(amr_wind::FieldState::NPH);
    icns().compute_diffusion_term(amr_wind::FieldState::Old);

    // *************************************************************************************
    // Evaluate right hand side and store in velocity
    // *************************************************************************************
    icns().compute_predictor_rhs(m_diff_type);

    if (m_verbose > 2) {
        PrintMaxVelLocations("after predictor rhs");
    }

    // *************************************************************************************
    // Solve diffusion equation for u* but using eta_old at old time
    // *************************************************************************************
    if (m_diff_type == DiffusionType::Crank_Nicolson ||
        m_diff_type == DiffusionType::Implicit) {
        amrex::Real dt_diff = (m_diff_type == DiffusionType::Implicit)
                                  ? m_time.delta_t()
                                  : 0.5_rt * m_time.delta_t();
        icns().solve(dt_diff);
    } else if (m_diff_type == DiffusionType::Explicit && m_use_godunov) {
        // explicit RK2
        auto diff_old = m_repo.create_scratch_field(
            AMREX_SPACEDIM, 0, amr_wind::FieldLoc::CELL);
        auto& diff_new =
            icns().fields().diff_term.state(amr_wind::FieldState::New);
        amr_wind::field_ops::copy(*diff_old, diff_new, 0, 0, AMREX_SPACEDIM, 0);
        icns().compute_diffusion_term(amr_wind::FieldState::New);
        amr_wind::field_ops::saxpy(
            diff_new, -1.0_rt, *diff_old, 0, 0, AMREX_SPACEDIM, 0);
        icns().improve_explicit_diffusion(m_time.delta_t());
    }
    icns().post_solve_actions();

    if (m_verbose > 2) {
        PrintMaxVelLocations("after diffusion solve");
    }

    // ************************************************************************************
    //
    // Project velocity field, update pressure
    //
    // ************************************************************************************
    ApplyProjection(
        (density_new).vec_const_ptrs(), new_time, m_time.delta_t(),
        incremental_projection);

    if (m_verbose > 2) {
        PrintMaxVelLocations("after nodal projection");
    }
    // ScratchField to store the old np1
    if (fixed_point_iteration > 0 && m_verbose > 0) {
        auto vel_np1_old = m_repo.create_scratch_field(
            "vel_np1_old", AMREX_SPACEDIM, 1, amr_wind::FieldLoc::CELL);

        auto vel_diff = m_repo.create_scratch_field(
            "vel_diff", AMREX_SPACEDIM, 1, amr_wind::FieldLoc::CELL);
        // lincomb to get old np1
        amr_wind::field_ops::lincomb(
            (*vel_np1_old), -1.0_rt,
            icns().fields().field.state(amr_wind::FieldState::Old), 0, 2,
            icns().fields().field.state(amr_wind::FieldState::NPH), 0, 0,
            icns().fields().field.num_comp(), 1);

        amr_wind::io::print_nonlinear_residual(m_sim, *vel_diff, *vel_np1_old);
        vel_np1_old.reset();
        vel_diff.reset();
    }
}
//
// Apply corrector:
//
//  Output variables from the predictor are labelled _pred
//
//  1. Use u = vel_pred to compute
//
//      conv_u  = - u grad u
//      conv_r  = - u grad rho
//      conv_t  = - u grad trac
//      eta     = viscosity
//      divtau  = div( eta ( (grad u) + (grad u)^T ) ) / rho
//
//      conv_u  = 0.5_rt (conv_u + conv_u_pred)
//      conv_r  = 0.5_rt (conv_r + conv_r_pred)
//      conv_t  = 0.5_rt (conv_t + conv_t_pred)
//      if (m_diff_type == DiffusionType::Explicit)
//         divtau  = divtau at new_time using (*) state
//      else
//         divtau  = 0.0_rt
//      eta     = eta at new_time
//
//     rhs = u + dt * ( conv + divtau )
//
//  2. Add explicit forcing term i.e. gravity + lagged pressure gradient
//
//      rhs += dt * ( g - grad(p + p0) / rho )
//
//      Note that in order to add the pressure gradient terms divided by rho,
//      we convert the velocity to momentum before adding and then convert them
//      back.
//
//  3. A. If (m_diff_type == DiffusionType::Implicit)
//        solve implicit diffusion equation for u*
//
//     ( 1 - dt / rho * div ( eta grad ) ) u* = u^n + dt * conv_u
//                                                  + dt * ( g - grad(p + p0) /
//                                                  rho )
//
//     B. If (m_diff_type == DiffusionType::Crank-Nicolson)
//        solve semi-implicit diffusion equation for u*
//
//     ( 1 - (dt/2) / rho * div ( eta grad ) ) u* = u^n + dt * conv_u + (dt/2) /
//     rho * div (eta_old grad) u^n
//                                                      + dt * ( g - grad(p +
//                                                      p0) / rho )
//
//  4. Apply projection
//
//     Add pressure gradient term back to u*:
//
//      u** = u* + dt * grad p / rho
//
//     Solve Poisson equation for phi:
//
//     div( grad(phi) / rho ) = div( u** )
//
//     Update pressure:
//
//     p = phi / dt
//
//     Update velocity, now divergence free
//
//     vel = u** - dt * grad p / rho
//

/** Corrector step for MOL scheme
 *
 *  <ol>
 *  <li> Solve transport equation for momentum and scalars
 *
 *  \f{align}
 *  \left[1 - \kappa \frac{\Delta t}{\rho} \nabla \cdot \left( \mu
 *  \nabla \right)\right] u^{*} &= u^n - \Delta t C_u + (1 - \kappa)
 * \frac{\Delta t}{\rho} \nabla \cdot \left( \mu \nabla\right) u^{n} +
 * \frac{\Delta t}{\rho} \left( S_u - \nabla(p + p_0)\right) \\ \f}
 *
 *  where
 *  \f{align}
 *  \kappa = \begin{cases}
 *  0 & \text{Explicit} \\
 *  0.5_rt & \text{Crank-Nicolson} \\
 *  1 & \text{Implicit}
 *  \end{cases}
 *  \f}
 *
 *  <li> \ref incflo::ApplyProjection "Apply projection"
 *  </ol>
 */
void incflo::ApplyCorrector()
{
    BL_PROFILE("amr-wind::incflo::ApplyCorrector");

    // We use the new time value for things computed on the "*" state
    amrex::Real new_time = m_time.new_time();

    if (m_verbose > 2) {
        PrintMaxValues("before corrector step");
    }

    amr_wind::io::print_mlmg_header("Corrector:");

    auto& density_new = density();
    const auto& density_old = density_new.state(amr_wind::FieldState::Old);
    auto& density_nph = density_new.state(amr_wind::FieldState::NPH);

    // Extrapolate and apply MAC projection for advection velocities
    icns().pre_advection_actions(amr_wind::FieldState::New);

    // *************************************************************************************
    // Compute the explicit "new" advective terms R_u^(n+1,*), R_r^(n+1,*) and
    // R_t^(n+1,*) We only reach the corrector if !m_use_godunov which means we
    // don't use the forces in constructing the advection term
    // *************************************************************************************
    for (auto& seqn : scalar_eqns()) {
        seqn->compute_advection_term(amr_wind::FieldState::New);
    }
    icns().compute_advection_term(amr_wind::FieldState::New);

    // *************************************************************************************
    // Compute viscosity / diffusive coefficients
    // *************************************************************************************
    m_sim.turbulence_model().update_turbulent_viscosity(
        amr_wind::FieldState::New, m_diff_type);
    icns().compute_mueff(amr_wind::FieldState::New);
    for (auto& eqns : scalar_eqns()) {
        eqns->compute_mueff(amr_wind::FieldState::New);
    }

    // Here we create divtau of the (n+1,*) state that was computed in the
    // predictor;
    //      we use this laps only if DiffusionType::Explicit
    if (m_diff_type == DiffusionType::Explicit) {
        icns().compute_diffusion_term(amr_wind::FieldState::New);

        for (auto& eqns : scalar_eqns()) {
            eqns->compute_diffusion_term(amr_wind::FieldState::New);
        }
    }

    // *************************************************************************************
    // Update density first
    // *************************************************************************************
    if (m_sim.pde_manager().constant_density()) {
        amr_wind::field_ops::copy(density_nph, density_old, 0, 0, 1, 1);
    }

    // TODO: This sub-section has not been adjusted for mesh mapping - adjust in
    // corrector too Perform scalar update one at a time. This is to allow an
    // updated density at `n+1/2` to be computed before other scalars use it
    // when computing their source terms.
    for (auto& eqn : scalar_eqns()) {
        // Compute (recompute for Godunov) the scalar forcing terms
        // Note this is (rho * scalar) and not just scalar
        eqn->compute_source_term(amr_wind::FieldState::New);

        // Update (note that dtdt already has rho in it)
        // (rho trac)^new = (rho trac)^old + dt * (
        //                   div(rho trac u) + div (mu grad trac) + rho * f_t
        eqn->compute_corrector_rhs(m_diff_type);

        auto& field = eqn->fields().field;
        if (m_diff_type != DiffusionType::Explicit) {
            amrex::Real dt_diff = (m_diff_type == DiffusionType::Implicit)
                                      ? m_time.delta_t()
                                      : 0.5_rt * m_time.delta_t();

            // Solve diffusion eqn. and update of the scalar field
            eqn->solve(dt_diff);
        }
        eqn->post_solve_actions();

        // Update scalar at n+1/2
        amr_wind::field_ops::lincomb(
            field.state(amr_wind::FieldState::NPH), 0.5_rt,
            field.state(amr_wind::FieldState::Old), 0, 0.5_rt, field, 0, 0,
            field.num_comp(), 1);
    }

    // *************************************************************************************
    // Define the forcing terms to use in the final update (using half-time
    // density)
    // *************************************************************************************
    icns().compute_source_term(amr_wind::FieldState::New);

    // *************************************************************************************
    // Evaluate right hand side and store in velocity
    // *************************************************************************************
    icns().compute_corrector_rhs(m_diff_type);

    // *************************************************************************************
    //
    // Solve diffusion equation for u* at t^{n+1} but using eta at predicted new
    // time
    //
    // *************************************************************************************

    if (m_diff_type == DiffusionType::Crank_Nicolson ||
        m_diff_type == DiffusionType::Implicit) {
        amrex::Real dt_diff = (m_diff_type == DiffusionType::Implicit)
                                  ? m_time.delta_t()
                                  : 0.5_rt * m_time.delta_t();
        icns().solve(dt_diff);
    }
    icns().post_solve_actions();

    // *************************************************************************************
    // Project velocity field, update pressure
    // *************************************************************************************
    bool incremental = false;
    ApplyProjection(
        (density_new).vec_const_ptrs(), new_time, m_time.delta_t(),
        incremental);
}

void incflo::prescribe_advance()
{
    BL_PROFILE("amr-wind::incflo::prescribe_advance");

    prepare_time_step();

    ApplyPrescribeStep();
}

void incflo::ApplyPrescribeStep()
{
    BL_PROFILE("amr-wind::incflo::ApplyPrescribeStep");
    // The intent of this function is to see the effect of a prescribed
    // advection velocity:
    // - No source terms or viscous terms are used for icns
    // - The MAC velocity is prescribed before this routine, not calculated here
    // - The nodal projection is omitted

    if (m_verbose > 2) {
        PrintMaxValues("before prescribe step");
    }

    auto& density_new = density();
    const auto& density_old = density_new.state(amr_wind::FieldState::Old);
    auto& density_nph = density_new.state(amr_wind::FieldState::NPH);

    // Compute diffusive and source terms for scalars
    m_sim.turbulence_model().update_turbulent_viscosity(
        amr_wind::FieldState::Old, m_diff_type);
    for (auto& eqns : scalar_eqns()) {
        eqns->compute_mueff(amr_wind::FieldState::Old);
    }
    if (m_use_godunov) {
        for (auto& seqn : scalar_eqns()) {
            seqn->compute_source_term(amr_wind::FieldState::Old);
        }
    }

    if (need_divtau()) {
        // Compute explicit diffusive terms for scalars
        for (auto& eqn : scalar_eqns()) {
            auto& field = eqn->fields().field;
            // Reuse existing buffer to avoid creating new multifabs
            amr_wind::field_ops::copy(
                field, field.state(amr_wind::FieldState::Old), 0, 0,
                field.num_comp(), 1);

            eqn->compute_diffusion_term(amr_wind::FieldState::Old);

            if (m_use_godunov) {
                amr_wind::field_ops::add(
                    eqn->fields().src_term, eqn->fields().diff_term, 0, 0,
                    field.num_comp(), 0);
            }
        }
    }

    if (m_use_godunov) {
        const int nghost_force = 1;
        amrex::IntVect ng(nghost_force);
        for (auto& eqn : scalar_eqns()) {
            eqn->fields().src_term.fillpatch(m_time.current_time(), ng);
        }
    }

    // For scalars only first
    // *************************************************************************************
    // if ( m_use_godunov) Compute the explicit advective terms
    //                     R_u^(n+1/2), R_s^(n+1/2) and R_t^(n+1/2)
    // if (!m_use_godunov) Compute the explicit advective terms
    //                     R_u^n      , R_s^n       and R_t^n
    // *************************************************************************************
    for (auto& seqn : scalar_eqns()) {
        seqn->compute_advection_term(amr_wind::FieldState::Old);
    }

    // *************************************************************************************
    // Update density first
    // *************************************************************************************
    if (m_sim.pde_manager().constant_density()) {
        amr_wind::field_ops::copy(density_nph, density_old, 0, 0, 1, 1);
    }

    for (auto& eqn : scalar_eqns()) {
        // Compute (recompute for Godunov) the scalar forcing terms
        eqn->compute_source_term(amr_wind::FieldState::NPH);

        // Update the scalar (if explicit), or the RHS for implicit/CN
        eqn->compute_predictor_rhs(m_diff_type);

        auto& field = eqn->fields().field;
        if (m_diff_type != DiffusionType::Explicit) {
            amrex::Real dt_diff = (m_diff_type == DiffusionType::Implicit)
                                      ? m_time.delta_t()
                                      : 0.5_rt * m_time.delta_t();

            // Solve diffusion eqn. and update of the scalar field
            eqn->solve(dt_diff);
        }
        // Post-processing actions after a PDE solve
        eqn->post_solve_actions();

        // Update scalar at n+1/2
        amr_wind::field_ops::lincomb(
            field.state(amr_wind::FieldState::NPH), 0.5_rt,
            field.state(amr_wind::FieldState::Old), 0, 0.5_rt, field, 0, 0,
            field.num_comp(), 1);
    }

    // With scalars computed, compute advection of momentum
    icns().compute_advection_term(amr_wind::FieldState::Old);

    // Evaluate right hand side and store in velocity
    // Explicit is used because viscous icns terms are supposed to be ignored
    icns().compute_predictor_rhs(DiffusionType::Explicit);

    icns().post_solve_actions();
}
