#include "src/equation_systems/icns/source_terms/PartialDragForcing.H"
#include "src/equation_systems/vof/volume_fractions.H"
#include "src/utilities/IOManager.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "src/wind_energy/ABL.H"
#include "src/physics/PartialTerrainDrag.H"
#include "src/utilities/linear_interpolation.H"
#include "src/utilities/constants.H"
#include "AMReX_REAL.H"
#include <fstream>

using namespace amrex::literals;

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::GpuArray<amrex::Real, 2>
compute_target_wind(
    const amrex::Real ux,
    const amrex::Real uy,
    const amrex::Real dx,
    const amrex::Real z0,
    const amrex::Real kappa)
{
    const amrex::Real wspd = std::sqrt(ux * ux + uy * uy);
    const amrex::Real ustar = wspd * kappa / std::log(1.5_rt * dx / z0);
    const amrex::Real wspd_target = ustar / kappa * std::log(0.5_rt * dx / z0);
    const amrex::Real wspd_target_x =
        wspd_target * ux /
        (kynema_sgf::constants::EPS + std::sqrt(ux * ux + uy * uy));
    const amrex::Real wspd_target_y =
        wspd_target * uy /
        (kynema_sgf::constants::EPS + std::sqrt(ux * ux + uy * uy));
    return {wspd_target_x, wspd_target_y};
}

} // namespace

namespace kynema_sgf::pde::icns {

PartialDragForcing::PartialDragForcing(const CFDSim& sim)
    : m_time(sim.time())
    , m_sim(sim)
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
{
    amrex::ParmParse pp("PartialDragForcing");
    pp.query("drag_coefficient", m_drag_coefficient);
    pp.query("bc_forcing_time_factor", m_bc_forcing_time_factor);
    pp.query("drag_forcing", m_drag_forcing);
    amrex::ParmParse pp_abl("ABL");
    pp_abl.query("wall_het_model", m_wall_het_model);
    pp_abl.query("monin_obukhov_length", m_monin_obukhov_length);
    pp_abl.query("kappa", m_kappa);
    pp_abl.query("mo_gamma_m", m_gamma_m);
    pp_abl.query("mo_beta_m", m_beta_m);
}

PartialDragForcing::~PartialDragForcing() = default;

void PartialDragForcing::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    auto const& src_arrs = src_term.arrays();
    auto const& vel =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_arrays();

    const int is_partial =
        this->m_sim.repo().field_exists("partial_terrain_blank") ? 1 : 0;
    if (is_partial == 0) {
        amrex::Abort(
            "Need partial terrain blanking variable to use this source term");
    }
    auto const& blank = this->m_sim.repo()
                            .get_field("partial_terrain_blank")(lev)
                            .const_arrays();

    const int has_terrainz0 =
        this->m_sim.repo().field_exists("partial_terrainz0") ? 1 : 0;
    auto const& terrainz0_arrs = has_terrainz0 != 0
                                     ? this->m_sim.repo()
                                           .get_field("partial_terrainz0")(lev)
                                           .const_arrays()
                                     : amrex::MultiArray4<amrex::Real const>();

    const auto& geom = m_mesh.Geom(lev);
    const auto& dx = geom.CellSizeArray();
    const amrex::Real drag_coefficient = m_drag_coefficient;
    const auto& dt = m_time.delta_t();
    const amrex::Real time_factor = m_bc_forcing_time_factor * dt;
    const amrex::Real Cd = drag_coefficient / dx[2];
    const amrex::Real kappa = m_kappa;
    const amrex::Real cd_max = 1000.0_rt;
    const int drag_forcing = m_drag_forcing;

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int /*n*/) {
            const amrex::Real ux1 = vel[nbx](i, j, k, 0);
            const amrex::Real uy1 = vel[nbx](i, j, k, 1);
            const amrex::Real uz1 = vel[nbx](i, j, k, 2);
            const amrex::Real m =
                std::sqrt((ux1 * ux1) + (uy1 * uy1) + (uz1 * uz1));
            amrex::Real bc_forcing_x = 0.0_rt;
            amrex::Real bc_forcing_y = 0.0_rt;
            amrex::Real bc_forcing_z = 0.0_rt;

            auto fluid_frac = [=](int ii, int jj, int kk) {
                return 1.0_rt - blank[nbx](ii, jj, kk, 0);
            };

            const amrex::Real z0 =
                has_terrainz0 != 0
                    ? amrex::max<amrex::Real>(
                          terrainz0_arrs[nbx](i, j, k), 1.0e-4_rt)
                    : 0.1_rt;

            //! West
            amrex::GpuArray<amrex::Real, 2> tmp_wind_target =
                compute_target_wind(
                    vel[nbx](i - 1, j, k, 2), vel[nbx](i - 1, j, k, 1), dx[0],
                    z0, kappa);
            bc_forcing_z += -(tmp_wind_target[0] - uz1) / time_factor *
                            fluid_frac(i + 1, j, k);
            bc_forcing_y += -(tmp_wind_target[1] - uy1) / time_factor *
                            fluid_frac(i + 1, j, k);
            //! East
            tmp_wind_target = compute_target_wind(
                vel[nbx](i + 1, j, k, 2), vel[nbx](i + 1, j, k, 1), dx[0], z0,
                kappa);
            bc_forcing_z += -(tmp_wind_target[0] - uz1) / time_factor *
                            fluid_frac(i - 1, j, k);
            bc_forcing_y += -(tmp_wind_target[1] - uy1) / time_factor *
                            fluid_frac(i - 1, j, k);
            //! South
            tmp_wind_target = compute_target_wind(
                vel[nbx](i, j - 1, k, 2), vel[nbx](i, j - 1, k, 0), dx[1], z0,
                kappa);
            bc_forcing_z += -(tmp_wind_target[0] - uz1) / time_factor *
                            fluid_frac(i, j + 1, k);
            bc_forcing_x += -(tmp_wind_target[1] - ux1) / time_factor *
                            fluid_frac(i, j + 1, k);
            //! North
            tmp_wind_target = compute_target_wind(
                vel[nbx](i, j + 1, k, 2), vel[nbx](i, j + 1, k, 0), dx[1], z0,
                kappa);
            bc_forcing_z += -(tmp_wind_target[0] - uz1) / time_factor *
                            fluid_frac(i, j - 1, k);
            bc_forcing_x += -(tmp_wind_target[1] - ux1) / time_factor *
                            fluid_frac(i, j - 1, k);
            //! Top
            tmp_wind_target = compute_target_wind(
                vel[nbx](i, j, k + 1, 0), vel[nbx](i, j, k + 1, 1), dx[2], z0,
                kappa);
            bc_forcing_x += -(tmp_wind_target[0] - ux1) / time_factor *
                            fluid_frac(i, j, k - 1);
            bc_forcing_y += -(tmp_wind_target[1] - uy1) / time_factor *
                            fluid_frac(i, j, k - 1);
            //! Bottom
            tmp_wind_target = compute_target_wind(
                vel[nbx](i, j, k - 1, 0), vel[nbx](i, j, k - 1, 1), dx[2], z0,
                kappa);
            bc_forcing_x += -(tmp_wind_target[0] - ux1) / time_factor *
                            fluid_frac(i, j, k + 1);
            bc_forcing_y += -(tmp_wind_target[1] - uy1) / time_factor *
                            fluid_frac(i, j, k + 1);

            const amrex::Real sum_blank_x =
                fluid_frac(i, j - 1, k) + fluid_frac(i, j + 1, k) +
                fluid_frac(i, j, k - 1) + fluid_frac(i, j, k + 1);
            bc_forcing_x /= (sum_blank_x + kynema_sgf::constants::EPS);
            const amrex::Real sum_blank_y =
                fluid_frac(i - 1, j, k) + fluid_frac(i + 1, j, k) +
                fluid_frac(i, j, k - 1) + fluid_frac(i, j, k + 1);
            bc_forcing_y /= (sum_blank_y + kynema_sgf::constants::EPS);
            const amrex::Real sum_blank_z =
                fluid_frac(i - 1, j, k) + fluid_frac(i + 1, j, k) +
                fluid_frac(i, j - 1, k) + fluid_frac(i, j + 1, k);
            bc_forcing_z /= (sum_blank_z + kynema_sgf::constants::EPS);

            // Target velocity intended for within terrain
            amrex::Real target_u = 0.;
            amrex::Real target_v = 0.;
            amrex::Real target_w = 0.;
            const amrex::Real CdM =
                std::min(Cd / (m + kynema_sgf::constants::EPS), cd_max / dx[2]);

            src_arrs[nbx](i, j, k, 0) -=
                (CdM * m * (ux1 - target_u) * fluid_frac(i, j, k) +
                 drag_forcing * bc_forcing_x * (1.0_rt - fluid_frac(i, j, k)));
            src_arrs[nbx](i, j, k, 1) -=
                (CdM * m * (uy1 - target_v) * fluid_frac(i, j, k) +
                 drag_forcing * bc_forcing_y * (1.0_rt - fluid_frac(i, j, k)));
            src_arrs[nbx](i, j, k, 2) -=
                (CdM * m * (uz1 - target_w) * fluid_frac(i, j, k) +
                 drag_forcing * bc_forcing_z * (1.0_rt - fluid_frac(i, j, k)));
        });
}

} // namespace kynema_sgf::pde::icns
