#include "amr-wind/wind_energy/ABL.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

ABLAnelastic::ABLAnelastic(CFDSim& sim) : m_sim(sim)
{
    {
        amrex::ParmParse pp("ABL");
        pp.query("anelastic", m_is_anelastic);
        pp.query("bottom_reference_pressure", m_bottom_reference_pressure);
    }
    std::string godunov_type;
    int conserv = 1;
    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
        pp.query("density", m_reference_density_constant);
        pp.query("godunov_type", godunov_type);
        pp.query("icns_conserv", conserv);
    }
    if (m_is_anelastic) {
        {
            // ensure use of perturbational pressure form
            amrex::ParmParse pp("ICNS");
            pp.add("use_perturb_pressure", true);
        }
        if (conserv != 1) {
            amrex::Abort("ABLAnelastic is not supported for icns_conserv != 1");
        }
        if (amrex::toLower(godunov_type) == "bds") {
            amrex::Abort("ABLAnelastic is not supported by BDS");
        }
        const auto& density = m_sim.repo().get_field("density");
        auto& ref_density = m_sim.repo().declare_field(
            "reference_density", 1, density.num_grow()[0], 1);
        ref_density.set_default_fillpatch_bc(m_sim.time());
        m_sim.repo().declare_nd_field("reference_pressure", 1, 0, 1);
        auto& ref_theta =
            m_sim.repo().declare_field("reference_temperature", 1, 1, 1);
        ref_theta.set_default_fillpatch_bc(m_sim.time());
    }
}

void ABLAnelastic::post_init_actions()
{
    if (!m_is_anelastic) {
        return;
    }
    initialize_data();
}

void ABLAnelastic::post_regrid_actions()
{
    if (!m_is_anelastic) {
        return;
    }
    initialize_data();
}

void ABLAnelastic::initialize_data()
{

    m_density.resize(m_axis, m_sim.mesh().Geom());
    m_pressure.resize(m_axis, m_sim.mesh().Geom());
    m_theta.resize(m_axis, m_sim.mesh().Geom());

    AMREX_ALWAYS_ASSERT(m_sim.repo().num_active_levels() == m_density.size());
    AMREX_ALWAYS_ASSERT(m_sim.repo().num_active_levels() == m_pressure.size());
    AMREX_ALWAYS_ASSERT(m_sim.repo().num_active_levels() == m_theta.size());

    initialize_isentropic_hse();

    auto& density0_field = m_sim.repo().get_field("reference_density");
    auto& p0 = m_sim.repo().get_field("reference_pressure");
    auto& temp0 = m_sim.repo().get_field("reference_temperature");
    m_density.copy_to_field(density0_field);
    m_pressure.copy_to_field(p0);
    m_theta.copy_to_field(temp0);
    density0_field.fillpatch(m_sim.time().current_time());
    temp0.fillpatch(m_sim.time().current_time());
    auto& density_field = m_sim.repo().get_field("density");
    field_ops::copy(
        density_field, density0_field, 0, 0, density0_field.num_comp(),
        density0_field.num_grow());
}

void ABLAnelastic::initialize_isentropic_hse()
{
    const int max_iterations = 10;
    const auto ref_theta = m_sim.transport_model().reference_temperature();
    const auto eos = amr_wind::eos::GammaLaw(m_bottom_reference_pressure);

    for (int lev = 0; lev < m_density.size(); lev++) {
        auto& dens = m_density.host_data(lev);
        auto& pres = m_pressure.host_data(lev);
        auto& theta = m_theta.host_data(lev);
        theta.assign(theta.size(), ref_theta);

        const auto dx = m_sim.mesh().Geom(lev).CellSize(m_axis);
        const amrex::Real half_dx = 0.5_rt * dx;

        // Initial guess
        dens[0] = m_reference_density_constant;
        pres[0] = m_bottom_reference_pressure +
                  (half_dx * dens[0] * m_gravity[m_axis]);

        // We do a Newton iteration to satisfy the EOS & HSE (with constant
        // theta) at the surface
        bool converged_hse = false;
        amrex::Real p_hse = 0.0_rt;
        amrex::Real p_eos = 0.0_rt;

        for (int iter = 0; (iter < max_iterations) && (!converged_hse);
             iter++) {
            p_hse = m_bottom_reference_pressure +
                    (half_dx * dens[0] * m_gravity[m_axis]);
            p_eos = eos.p_rth(dens[0], ref_theta);

            const amrex::Real p_diff = p_hse - p_eos;
            const amrex::Real dpdr = eos.dp_constanttheta(dens[0], ref_theta);
            const amrex::Real ddens =
                p_diff / (dpdr - half_dx * m_gravity[m_axis]);

            dens[0] = dens[0] + ddens;
            pres[0] = eos.p_rth(dens[0], ref_theta);

            if (std::abs(ddens) < constants::LOOSE_TOL) {
                converged_hse = true;
                break;
            }
        }
        if (!converged_hse) {
            amrex::Abort("Initializing with hydrostatic equilibrium failed");
        }

        // To get values at k > 0 we do a Newton iteration to satisfy the EOS
        // (with constant theta)
        for (int k = 1; k < dens.size(); k++) {
            converged_hse = false;

            dens[k] = dens[k - 1];
            for (int iter = 0; (iter < max_iterations) && (!converged_hse);
                 iter++) {
                const amrex::Real dens_avg = 0.5_rt * (dens[k - 1] + dens[k]);
                p_hse = pres[k - 1] + (dx * dens_avg * m_gravity[m_axis]);
                p_eos = eos.p_rth(dens[k], ref_theta);

                const amrex::Real p_diff = p_hse - p_eos;
                const amrex::Real dpdr =
                    eos.dp_constanttheta(dens[k], ref_theta);
                const amrex::Real ddens =
                    p_diff / (dpdr - dx * m_gravity[m_axis]);

                dens[k] = dens[k] + ddens;
                pres[k] = eos.p_rth(dens[k], ref_theta);

                if (std::abs(ddens) < constants::LOOSE_TOL * dens[k - 1]) {
                    converged_hse = true;
                    break;
                }
            }
            if (!converged_hse) {
                amrex::Abort(
                    "Initializing with hydrostatic equilibrium failed");
            }
        }
    }
    m_density.copy_host_to_device();
    m_pressure.copy_host_to_device();
    m_theta.copy_host_to_device();
}

} // namespace amr_wind
