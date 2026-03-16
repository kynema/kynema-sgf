#include "amr-wind/overset/OversetOps.H"
#include "amr-wind/overset/overset_ops_routines.H"
#include "amr-wind/core/field_ops.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFabUtil.H"
#include "amr-wind/core/MLMGOptions.H"
#include "amr-wind/projection/nodal_projection_ops.H"
#include <hydro_NodalProjector.H>
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/wind_energy/ABLBoundaryPlane.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

void OversetOps::initialize(CFDSim& sim)
{
    m_sim_ptr = &sim;
    // Queries for reinitialization options
    amrex::ParmParse pp("Overset");

    // Queries for coupling options
    pp.query("replace_gradp_postsolve", m_replace_gradp_postsolve);
    pp.query("verbose", m_verbose);

    m_vof_exists = m_sim_ptr->repo().field_exists("vof");
    if (m_vof_exists) {
        m_mphase = &m_sim_ptr->physics_manager().get<MultiPhase>();

        // Check combination of pressure options
        if (m_mphase->perturb_pressure() &&
            !m_mphase->reconstruct_true_pressure()) {
            amrex::Abort(
                "OversetOps: perturbational pressure is turned on, but true "
                "pressure reconstruction is turned off. This approach will be "
                "incorrect when coupling with Nalu-Wind.");
        }
        if (!m_mphase->perturb_pressure() &&
            m_mphase->reconstruct_true_pressure()) {
            amrex::Print()
                << "WARNING (OversetOps): true pressure reconstruction is "
                   "turned on, but it will remain inactive because "
                   "perturbational pressure is turned off.\n";
        }
        if (m_replace_gradp_postsolve) {
            m_gp_copy = &m_sim_ptr->repo().declare_field("gp_copy", 3);
        }
    }

    parameter_output();
}

void OversetOps::pre_advance_work()
{
    // Update pressure gradient using updated overset pressure field
    update_gradp();

    if (m_vof_exists) {
        if (m_mphase->perturb_pressure()) {
            // Modify to be consistent with internal source terms
            form_perturb_pressure();
        }
        // Update pressure gradient using sharpened pressure field
        update_gradp();
        // Calculate vof-dependent node mask
        const auto& iblank = m_sim_ptr->repo().get_int_field("iblank_node");
        const auto& vof = m_sim_ptr->repo().get_field("vof");
        auto& mask = m_sim_ptr->repo().get_int_field("mask_node");
        overset_ops::iblank_node_to_mask_vof(iblank, vof, mask);

        // If pressure gradient will be replaced, store current pressure
        // gradient
        if (m_replace_gradp_postsolve) {
            const auto& gp = m_sim_ptr->repo().get_field("gp");
            for (int lev = 0; lev < m_sim_ptr->repo().num_active_levels();
                 ++lev) {
                amrex::MultiFab::Copy(
                    (*m_gp_copy)(lev), gp(lev), 0, 0, gp(lev).nComp(),
                    (m_gp_copy)->num_grow());
            }
        }
    }

    // Pre advance work for plane was skipped for overset solver, do it here
    if (m_sim_ptr->physics_manager().contains("ABL")) {
        auto& abl = m_sim_ptr->physics_manager().get<ABL>();
        abl.bndry_plane().pre_advance_work();
    }
}

void OversetOps::update_gradp()
{
    BL_PROFILE("amr-wind::OversetOps::update_gradp");

    // Get relevant fields
    auto& grad_p = m_sim_ptr->repo().get_field("gp");
    auto& pressure = m_sim_ptr->repo().get_field("p");
    auto& velocity = m_sim_ptr->repo().get_field("velocity");

    // Set up projection object
    std::unique_ptr<Hydro::NodalProjector> nodal_projector;
    // Boundary conditions from field, periodicity
    auto bclo = nodal_projection::get_projection_bc(
        amrex::Orientation::low, pressure,
        m_sim_ptr->mesh().Geom(0).isPeriodic());
    auto bchi = nodal_projection::get_projection_bc(
        amrex::Orientation::high, pressure,
        m_sim_ptr->mesh().Geom(0).isPeriodic());
    // Velocity multifab is needed for proper initialization, but only the size
    // matters for the purpose of calculating gradp, the values do not matter
    int finest_level = m_sim_ptr->mesh().finestLevel();
    amrex::Vector<amrex::MultiFab*> vel;
    for (int lev = 0; lev <= finest_level; ++lev) {
        vel.push_back(&(velocity(lev)));
    }
    amr_wind::MLMGOptions options("nodal_proj");
    // Create nodal projector with unity scaling factor for simplicity
    nodal_projector = std::make_unique<Hydro::NodalProjector>(
        vel, 1.0_rt, m_sim_ptr->mesh().Geom(0, finest_level), options.lpinfo());
    // Set MLMG and NodalProjector options
    options(*nodal_projector);
    nodal_projector->setDomainBC(bclo, bchi);

    // Recalculate gradphi with fluxes
    auto gradphi = nodal_projector->calcGradPhi(pressure.vec_ptrs());

    // Transfer pressure gradient to gp field
    for (int lev = 0; lev <= finest_level; lev++) {
        const auto& gp_lev_arrs = grad_p(lev).arrays();
        const auto& gp_proj_arrs = gradphi[lev]->const_arrays();
        amrex::ParallelFor(
            grad_p(lev), amrex::IntVect(0), AMREX_SPACEDIM,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
                gp_lev_arrs[nbx](i, j, k, n) = gp_proj_arrs[nbx](i, j, k, n);
            });
    }
    amrex::Gpu::streamSynchronize();

    // Averaging down here would be unnecessary; it is built into calcGradPhi
}

void OversetOps::post_advance_work()
{
    // Replace and reapply pressure gradient if requested
    if (m_vof_exists && m_replace_gradp_postsolve) {
        replace_masked_gradp();
    }
}

/* ----------------------------------------------- */
/* PUBLIC FUNCTIONS ABOVE, PRIVATE FUNCTIONS BELOW */
/* ----------------------------------------------- */

void OversetOps::parameter_output() const
{
    // Print the details
    if (m_verbose > 0 && m_vof_exists) {
        // Important parameters
        amrex::Print() << "\nOverset Coupling Parameters: \n"
                       << "---- Replace overset pres grad: "
                       << m_replace_gradp_postsolve << '\n';
        amrex::Print() << "---- Perturbational pressure  : "
                       << m_mphase->perturb_pressure() << '\n'
                       << "---- Reconstruct true pressure: "
                       << m_mphase->reconstruct_true_pressure() << std::endl;
        amrex::Print() << std::endl;
    }
}

void OversetOps::form_perturb_pressure()
{
    auto& pressure = m_sim_ptr->repo().get_field("p");
    const auto& p0 = m_sim_ptr->repo().get_field("reference_pressure");
    for (int lev = 0; lev < m_sim_ptr->repo().num_active_levels(); lev++) {
        amrex::MultiFab::Subtract(
            pressure(lev), p0(lev), 0, 0, 1, pressure.num_grow()[0]);
    }
}

void OversetOps::replace_masked_gradp()
{
    const auto& repo = m_sim_ptr->repo();
    const auto nlevels = repo.num_active_levels();

    // Get timestep
    const amrex::Real dt = m_sim_ptr->time().delta_t();

    // Get blanking for cells
    const auto& iblank_cell = repo.get_int_field("iblank_cell");

    // Get fields that will be modified or used
    auto& vel = repo.get_field("velocity");
    auto& rho = repo.get_field("density");
    auto& gp = repo.get_field("gp");

    // For iblanked cells, replace gp with original gp, to get original vel
    for (int lev = 0; lev < nlevels; ++lev) {
        // Remove pressure gradient term
        overset_ops::apply_pressure_gradient(vel(lev), rho(lev), gp(lev), -dt);
        // Modify pressure gradient
        overset_ops::replace_gradp(
            gp(lev), (*m_gp_copy)(lev), iblank_cell(lev));
        // Reapply pressure gradient term
        overset_ops::apply_pressure_gradient(vel(lev), rho(lev), gp(lev), dt);
    }
    amrex::Gpu::streamSynchronize();
}

} // namespace amr_wind
