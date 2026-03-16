#include "amr-wind/equation_systems/PDEBase.H"
#include "amr-wind/equation_systems/SchemeTraits.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/equation_systems/PDEHelpers.H"
#include "amr-wind/incflo_enums.H"
#include "amr-wind/core/field_ops.H"
#include "amr-wind/utilities/constants.H"
#include "amr-wind/utilities/diagnostics.H"

#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::pde {

PDEFields::PDEFields(FieldRepo& repo_in, const std::string& var_name)
    : repo(repo_in)
    , field(repo.get_field(var_name))
    , mueff(repo.get_field(pde_impl::mueff_name(var_name)))
    , src_term(repo.get_field(pde_impl::src_term_name(var_name)))
    , diff_term(repo.get_field(pde_impl::diff_term_name(var_name)))
    , conv_term(repo.get_field(pde_impl::conv_term_name(var_name)))
{}

PDEMgr::PDEMgr(CFDSim& sim) : m_sim(sim)
{
    amrex::ParmParse pp("incflo");
    pp.query("use_godunov", m_use_godunov);
    pp.query("constant_density", m_constant_density);

    bool is_anelastic = false;
    {
        amrex::ParmParse pp_abl("ABL");
        pp_abl.query("anelastic", is_anelastic);
    }

    if ((is_anelastic) && (m_constant_density)) {
        amrex::Print()
            << "WARNING: Anelastic implies variable density. Setting "
               "contant_density to false. Add incflo.constant_density=false to "
               "your input file to remove this warning."
            << '\n';
        m_constant_density = false;
    }

    m_scheme =
        m_use_godunov ? fvm::Godunov::scheme_name() : fvm::MOL::scheme_name();
}

PDEBase& PDEMgr::register_icns()
{
    const std::string name = "ICNS-" + m_scheme;

    m_icns = amr_wind::pde::PDEBase::create(name, m_sim);
    return *m_icns;
}

PDEBase& PDEMgr::register_transport_pde(const std::string& pde_name)
{
    const std::string name = pde_name + "-" + m_scheme;

    if (contains(name)) {
        amrex::Print() << "WARNING: Multiple requests to register PDE: "
                       << pde_name << '\n';
        return operator()(name);
    }

    return create(name, m_sim);
}

bool PDEMgr::has_pde(const std::string& pde_name) const
{
    const std::string name = pde_name + "-" + m_scheme;
    return contains(name);
}

int PDEMgr::num_ghost_state() const
{
    return m_use_godunov ? fvm::Godunov::nghost_state : fvm::MOL::nghost_state;
}

void PDEMgr::advance_states()
{
    if (m_constant_density) {
        m_sim.repo().get_field("density").advance_states();
    }

    icns().fields().field.advance_states();
    for (auto& eqn : scalar_eqns()) {
        eqn->fields().field.advance_states();
    }
}

void PDEMgr::prepare_boundaries()
{
    // If state variables exist at NPH, fill their boundary cells
    const auto nph_time =
        0.5_rt * (m_sim.time().current_time() + m_sim.time().new_time());
    if (m_constant_density &&
        m_sim.repo().field_exists("density", FieldState::NPH)) {
        auto& nph_field =
            m_sim.repo().get_field("density").state(FieldState::NPH);
        auto& new_field =
            m_sim.repo().get_field("density").state(FieldState::New);
        // Need to check if this step is necessary
        // Guessing that for no-op dirichlet bcs, the new needs to be copied to
        // NPH
        field_ops::copy(
            nph_field, new_field, 0, 0, new_field.num_comp(),
            new_field.num_grow());
    }

    if (m_sim.repo().field_exists(
            icns().fields().field.base_name(), FieldState::NPH)) {
        auto& nph_field = icns().fields().field.state(FieldState::NPH);
        auto& new_field = icns().fields().field.state(FieldState::New);
        field_ops::copy(
            nph_field, new_field, 0, 0, new_field.num_comp(),
            new_field.num_grow());
        nph_field.fillphysbc(nph_time);
    }
    for (auto& eqn : scalar_eqns()) {
        eqn->fields().field.advance_states();
        if (m_sim.repo().field_exists(
                eqn->fields().field.base_name(), FieldState::NPH)) {
            auto& nph_field = eqn->fields().field.state(FieldState::NPH);
            auto& new_field = eqn->fields().field.state(FieldState::New);
            field_ops::copy(
                nph_field, new_field, 0, 0, new_field.num_comp(),
                new_field.num_grow());
            nph_field.fillphysbc(nph_time);
        }
    }

    // Fill mac velocities (ghost cells) using velocity BCs
    auto& u_mac = m_sim.repo().get_field("u_mac");
    auto& v_mac = m_sim.repo().get_field("v_mac");
    auto& w_mac = m_sim.repo().get_field("w_mac");
    const amrex::IntVect zeros{0, 0, 0};
    if (u_mac.num_grow() > zeros) {
        amrex::Array<Field*, AMREX_SPACEDIM> mac_vel = {
            AMREX_D_DECL(&u_mac, &v_mac, &w_mac)};
        icns().fields().field.fillpatch_sibling_fields(
            nph_time, u_mac.num_grow(), mac_vel);
    }
}

void PDEMgr::density_check()
{
    std::string advice(
        "\nCheck the initial bulk density and the inflow BC density values in "
        "the input file. When not specified, the default density is 1.0_rt.");
    std::string advice2(
        "\nIf this simulation begins from a restart file, confirm that the "
        "previous density is compatible with the parameters in the input "
        "file.");
    std::string advice3(
        "\nIf specific scalar boundary conditions are specified, make sure "
        "that vof and density BCs are the same (if Neumann type) or compatible "
        "(if Dirichlet type).");
    if (m_sim.repo().field_exists("vof")) {
        amrex::Real rho_l_max{1.0_rt}, rho_l_min{1.0_rt};
        amrex::Real rho_g_max{1.0_rt}, rho_g_min{1.0_rt};
        const bool liquid_found = diagnostics::get_field_extrema(
            rho_l_max, rho_l_min, m_sim.repo().get_field("density"),
            m_sim.repo().get_field("vof"), 1.0_rt, 0, 1, 1);
        const bool gas_found = diagnostics::get_field_extrema(
            rho_g_max, rho_g_min, m_sim.repo().get_field("density"),
            m_sim.repo().get_field("vof"), 0.0_rt, 0, 1, 1);
        if (liquid_found &&
            std::abs(rho_l_max - rho_l_min) > constants::LOOSE_TOL) {
            amrex::Abort(
                "Density check failed. Liquid density maximum is too different "
                "from liquid density minimum.\n"
                "rho_l_max = " +
                std::to_string(rho_l_max) + ", rho_l_min = " +
                std::to_string(rho_l_min) + advice2 + advice3);
        }
        if (gas_found &&
            std::abs(rho_g_max - rho_g_min) > constants::LOOSE_TOL) {
            amrex::Abort(
                "Density check failed. Gas density maximum is too different "
                "from gas density minimum.\n"
                "rho_g_max = " +
                std::to_string(rho_g_max) + ", rho_g_min = " +
                std::to_string(rho_g_min) + advice + advice2 + advice3);
        }
        if (!liquid_found) {
            amrex::Print()
                << "WARNING: Simulation is multiphase, but no fully liquid "
                   "cells (vof = 1) were found during the initial density "
                   "check.\n";
        }
        if (!gas_found) {
            amrex::Print()
                << "WARNING: Simulation is multiphase, but no fully gas "
                   "cells (vof = 0) were found during the initial density "
                   "check.\n";
        }
    } else if (m_constant_density) {
        amrex::Real rho_max{1.0_rt}, rho_min{1.0_rt};
        diagnostics::get_field_extrema(
            rho_max, rho_min, m_sim.repo().get_field("density"), 0, 1, 1);
        if (std::abs(rho_max - rho_min) > constants::LOOSE_TOL) {
            amrex::Abort(
                "Density check failed. Density maximum is too different "
                "from minimum for a constant density simulation.\n"
                "rho_max = " +
                std::to_string(rho_max) +
                ", rho_min = " + std::to_string(rho_min) + advice + advice2);
        }
    }
}

void PDEMgr::fillpatch_state_fields(
    const amrex::Real time, const FieldState fstate)
{
    if (m_constant_density) {
        m_sim.repo().get_field("density").state(fstate).fillpatch(time);
    }

    icns().fields().field.state(fstate).fillpatch(time);
    for (auto& eqn : scalar_eqns()) {
        eqn->fields().field.state(fstate).fillpatch(time);
    }
}

} // namespace amr_wind::pde
