#ifndef VOF_ADVECTION_H
#define VOF_ADVECTION_H

#include "amr-wind/equation_systems/vof/vof.H"
#include "amr-wind/equation_systems/vof/vof_hybridsolver_ops.H"
#include "amr-wind/equation_systems/vof/SplitAdvection.H"
#include "amr-wind/core/SimTime.H"

namespace amr_wind::pde {

/** Godunov scheme for VOF
 *  \ingroup vof
 */
template <>
struct AdvectionOp<VOF, fvm::Godunov>
{
    AdvectionOp(
        CFDSim& sim,
        PDEFields& fields_in,
        bool /*unused*/,
        bool /*unused*/,
        bool /*unused*/,
        bool /*unused*/)
        : m_time(sim.time())
        , fields(fields_in)
        , u_mac(fields_in.repo.get_field("u_mac"))
        , v_mac(fields_in.repo.get_field("v_mac"))
        , w_mac(fields_in.repo.get_field("w_mac"))
    {
        amrex::ParmParse pp_multiphase("VOF");
        pp_multiphase.query("remove_debris", m_rm_debris);
        pp_multiphase.query("replace_masked", m_replace_mask);

        // Setup density factor arrays for multiplying velocity flux
        fields_in.repo.declare_face_normal_field(
            {"advalpha_x", "advalpha_y", "advalpha_z"}, 1,
            fvm::Godunov::nghost_mac, 1);
    }

    void preadvect(
        const FieldState /*unused*/,
        const amrex::Real /*unused*/,
        const amrex::Real /*unused*/)
    {}

    void operator()(const FieldState /*unused*/, const amrex::Real dt)
    {
        static_assert(
            VOF::ndim == 1, "Invalid number of components for scalar");

        auto& repo = fields.repo;
        const auto& geom = repo.mesh().Geom();
        const int nlevels = repo.num_active_levels();

        auto& aa_x = repo.get_field("advalpha_x");
        auto& aa_y = repo.get_field("advalpha_y");
        auto& aa_z = repo.get_field("advalpha_z");

        // Old and new states
        auto& dof_old = fields.field.state(amr_wind::FieldState::Old);
        auto& dof_new = fields.field;
        // Working state of vof is nph, to keep others untouched during step
        auto& dof_field = fields.field.state(amr_wind::FieldState::NPH);

        // Initialize as old state values
        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            amrex::MultiFab::Copy(
                dof_field(lev), dof_old(lev), 0, 0, dof_field.num_comp(),
                dof_field.num_grow());
        }

        //
        // Advect volume using Implicit Eulerian Sweeping method with PLIC
        // reconstruction
        //

        auto flux_x =
            repo.create_scratch_field(1, 0, amr_wind::FieldLoc::XFACE);
        auto flux_y =
            repo.create_scratch_field(1, 0, amr_wind::FieldLoc::YFACE);
        auto flux_z =
            repo.create_scratch_field(1, 0, amr_wind::FieldLoc::ZFACE);

        // Scratch field for fluxC
        auto fluxC = repo.create_scratch_field(1, 0, amr_wind::FieldLoc::CELL);

        // Define the sweep time
        isweep += 1;
        if (isweep > 3) {
            isweep = 1;
        }

        amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> fluxes(
            nlevels);
        amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> advas(
            nlevels);
        for (int lev = 0; lev < nlevels; ++lev) {
            fluxes[lev][0] = &(*flux_x)(lev);
            fluxes[lev][1] = &(*flux_y)(lev);
            fluxes[lev][2] = &(*flux_z)(lev);
            advas[lev][0] = &aa_x(lev);
            advas[lev][1] = &aa_y(lev);
            advas[lev][2] = &aa_z(lev);
        }

        // Split advection step 1, with cmask calculation
        multiphase::split_advection_step(
            isweep, 0, nlevels, dof_field, fluxes, (*fluxC), advas, u_mac,
            v_mac, w_mac, dof_field.bc_type(), geom, m_time.new_time(), dt,
            m_rm_debris);
        // (copy old boundaries to working state)
        // Split advection step 2
        multiphase::split_advection_step(
            isweep, 1, nlevels, dof_field, fluxes, (*fluxC), advas, u_mac,
            v_mac, w_mac, dof_field.bc_type(), geom, m_time.new_time(), dt,
            m_rm_debris);
        // (copy old boundaries to working state)
        // Split advection step 3
        multiphase::split_advection_step(
            isweep, 2, nlevels, dof_field, fluxes, (*fluxC), advas, u_mac,
            v_mac, w_mac, dof_field.bc_type(), geom, m_time.new_time(), dt,
            m_rm_debris);

        // Replace masked cells using overset
        if (repo.int_field_exists("iblank_cell") && m_replace_mask) {
            auto& f_iblank = repo.get_int_field("iblank_cell");
            amr_wind::multiphase::replace_masked_vof(
                nlevels, f_iblank, dof_field, dof_new);
        }

        // Copy working version of vof to new state
        for (int lev = 0; lev < repo.num_active_levels(); ++lev) {
            amrex::MultiFab::Copy(
                dof_new(lev), dof_field(lev), 0, 0, dof_field.num_comp(),
                dof_field.num_grow());
        }
    }

    const SimTime& m_time;
    PDEFields& fields;
    Field& u_mac;
    Field& v_mac;
    Field& w_mac;
    int isweep = 0;
    bool m_rm_debris{true};
    bool m_replace_mask{true};
    // Lagrangian transport is deprecated, only Eulerian is supported
};

} // namespace amr_wind::pde
#endif
