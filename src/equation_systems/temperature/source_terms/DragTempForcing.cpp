#include "src/equation_systems/temperature/source_terms/DragTempForcing.H"
#include "src/utilities/IOManager.H"
#include "src/wind_energy/MOData.H"
#include "src/physics/TerrainDrag.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::temperature {

DragTempForcing::DragTempForcing(const CFDSim& sim)
    : m_time(sim.time())
    , m_sim(sim)
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_temperature(sim.repo().get_field("temperature"))
{
    amrex::ParmParse pp("DragTempForcing");
    pp.query("drag_coefficient", m_drag_coefficient);
    pp.query("soil_temperature", m_soil_temperature);
    pp.query("bc_forcing_time_factor", m_forcing_time_factor);
    // On with TerrainDrag.wall_treatment = improved, and can be set on its
    // own
    m_blank_follow_fluid = terraindrag::improved_wall_treatment();
    pp.query("blank_follow_fluid", m_blank_follow_fluid);
    amrex::ParmParse pp_abl("ABL");
    pp_abl.query("wall_het_model", m_wall_het_model);
    pp_abl.query("monin_obukhov_length", m_monin_obukhov_length);
    pp_abl.query("kappa", m_kappa);
    pp_abl.query("mo_gamma_m", m_gamma_m);
    pp_abl.query("mo_beta_m", m_beta_m);
    pp_abl.query("mo_gamma_m", m_gamma_h);
    pp_abl.query("mo_beta_m", m_beta_h);

    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.queryarr("gravity", m_gravity);
}

DragTempForcing::~DragTempForcing() = default;

void DragTempForcing::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const auto& vel_arrs =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_arrays();
    const auto& temperature_arrs =
        m_temperature.state(field_impl::dof_state(fstate))(lev).const_arrays();
    const bool is_terrain =
        this->m_sim.repo().int_field_exists("terrain_blank");
    if (!is_terrain) {
        amrex::Abort("Need terrain blanking variable to use this source term");
    }
    const auto& blank_arrs =
        this->m_sim.repo().get_int_field("terrain_blank")(lev).const_arrays();
    const auto& drag_arrs =
        this->m_sim.repo().get_int_field("terrain_drag")(lev).const_arrays();
    const auto& terrainz0_arrs =
        this->m_sim.repo().get_field("terrainz0")(lev).const_arrays();
    const auto& geom = m_mesh.Geom(lev);
    const auto& dx = geom.CellSizeArray();
    const amrex::Real drag_coefficient = m_drag_coefficient / dx[2];
    const amrex::Real gravity_mod = std::abs(m_gravity[2]);
    const amrex::Real kappa = m_kappa;
    const amrex::Real z0_min = 1.0e-4_rt;
    const amrex::Real monin_obukhov_length = m_monin_obukhov_length;
    const auto& dt = m_time.delta_t();
    const amrex::Real psi_m =
        (m_wall_het_model == "mol")
            ? MOData::calc_psi_m(
                  1.5_rt * dx[2] / m_monin_obukhov_length, m_beta_m, m_gamma_m)
            : 0.0_rt;
    const amrex::Real psi_h_neighbour =
        (m_wall_het_model == "mol")
            ? MOData::calc_psi_h(
                  1.5_rt * dx[2] / m_monin_obukhov_length, m_beta_h, m_gamma_h)
            : 0.0_rt;
    const amrex::Real psi_h_cell =
        (m_wall_het_model == "mol")
            ? MOData::calc_psi_h(
                  0.5_rt * dx[2] / m_monin_obukhov_length, m_beta_h, m_gamma_h)
            : 0.0_rt;
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real cd_max = 10.0_rt;
    const amrex::Real T0 = m_soil_temperature;
    const amrex::Real time_factor = m_forcing_time_factor;

    auto const& src_arrs = src_term.arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), 1,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
            const auto& vel = vel_arrs[nbx];
            const auto& temperature = temperature_arrs[nbx];
            const auto& blank = blank_arrs[nbx];
            const auto& drag = drag_arrs[nbx];
            const auto& terrainz0 = terrainz0_arrs[nbx];

            const amrex::Real z0 =
                amrex::max<amrex::Real>(terrainz0(i, j, k), z0_min);
            const amrex::Real ux1 = vel(i, j, k, 0);
            const amrex::Real uy1 = vel(i, j, k, 1);
            const amrex::Real uz1 = vel(i, j, k, 2);
            const amrex::Real theta = temperature(i, j, k, 0);
            const amrex::Real theta2 = temperature(i, j, k + 1, 0);
            const amrex::Real wspd = std::sqrt((ux1 * ux1) + (uy1 * uy1));
            const amrex::Real ustar =
                wspd * kappa / (std::log(1.5_rt * dx[2] / z0) - psi_m);
            const amrex::Real thetastar =
                theta * ustar * ustar /
                (kappa * gravity_mod * monin_obukhov_length);
            const amrex::Real surf_temp =
                theta2 - (thetastar / kappa *
                          (std::log(1.5_rt * dx[2] / z0) - psi_h_neighbour));
            const amrex::Real tTarget =
                surf_temp + (thetastar / kappa *
                             (std::log(0.5_rt * dx[2] / z0) - psi_h_cell));
            const amrex::Real bc_forcing_t =
                -(tTarget - theta) / (time_factor * dt);
            const amrex::Real m =
                std::sqrt((ux1 * ux1) + (uy1 * uy1) + (uz1 * uz1));
            const amrex::Real Cd = amrex::min<amrex::Real>(
                drag_coefficient / (m + tiny), cd_max / dx[2]);
            src_arrs[nbx](i, j, k, 0) -=
                ((Cd * (theta - T0) * blank(i, j, k, 0)) +
                 (bc_forcing_t * drag(i, j, k)));
        });

    if (m_blank_follow_fluid) {
        blank_follow_fluid(lev, fstate, src_term);
    }
}

// DragTempForcing.blank_follow_fluid: operator() relaxes the blanked cells
//  toward the soil temperature T0 at the rate C_d. Held at T0 under air of
//  another temperature, the blanked column conducts heat through its top face
//  into the drag cell and the cells beside it, a surface flux on top of the
//  one the drag-cell forcing sets. Here the blanked cells relax instead toward
//  the temperature of the cell above them, so that no gradient forms across
//  the blanked face, with the relaxation integrated exactly over the step:
//      S = -C_e (theta_k - theta_{k+1}),   C_e = (1 - exp(-C_d dt)) / dt,
//      C_d = min(drag_coefficient / (dz |u|), 10 / dz).
//  The explicit rate C_d reaches 2 per step at dz = 5 m and dt = 1 s, and the
//  coupled blanked column then grows without bound. The soil term that
//  operator() applied is removed first; the drag cells are not changed.
void DragTempForcing::blank_follow_fluid(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const auto& vel_arrs =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_arrays();
    const auto& temperature_arrs =
        m_temperature.state(field_impl::dof_state(fstate))(lev).const_arrays();
    const auto& blank_arrs =
        this->m_sim.repo().get_int_field("terrain_blank")(lev).const_arrays();
    const auto& dx = m_mesh.Geom(lev).CellSizeArray();
    const amrex::Real drag_coefficient = m_drag_coefficient / dx[2];
    const amrex::Real cd_max = 10.0_rt;
    const amrex::Real T0 = m_soil_temperature;
    const amrex::Real dt = m_time.delta_t();
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    auto const& src_arrs = src_term.arrays();
    amrex::ParallelFor(
        src_term, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            if (blank_arrs[nbx](i, j, k) != 1) {
                return;
            }
            const auto& vel = vel_arrs[nbx];
            const auto& temperature = temperature_arrs[nbx];
            const amrex::Real m = std::sqrt(
                (vel(i, j, k, 0) * vel(i, j, k, 0)) +
                (vel(i, j, k, 1) * vel(i, j, k, 1)) +
                (vel(i, j, k, 2) * vel(i, j, k, 2)));
            const amrex::Real Cd = amrex::min<amrex::Real>(
                drag_coefficient / (m + tiny), cd_max / dx[2]);
            const amrex::Real C_follow = (1.0_rt - std::exp(-Cd * dt)) / dt;
            const amrex::Real theta = temperature(i, j, k);
            src_arrs[nbx](i, j, k) +=
                (Cd * (theta - T0)) -
                (C_follow * (theta - temperature(i, j, k + 1)));
        });
}

} // namespace kynema_sgf::pde::temperature
