#include "amr-wind/equation_systems/icns/source_terms/DragForcing.H"
#include "amr-wind/equation_systems/vof/volume_fractions.H"
#include "amr-wind/utilities/IOManager.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "amr-wind/wind_energy/ABL.H"
#include "amr-wind/physics/TerrainDrag.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "amr-wind/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace {
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real viscous_drag_calculations(
    amrex::Real& Dxz,
    amrex::Real& Dyz,
    const amrex::Real ux1r,
    const amrex::Real uy1r,
    const amrex::Real ux2r,
    const amrex::Real uy2r,
    const amrex::Real z0,
    const amrex::Real dz,
    const amrex::Real kappa,
    const amrex::Real non_neutral_neighbour)
{
    const amrex::Real m2 = std::sqrt((ux2r * ux2r) + (uy2r * uy2r));
    const amrex::Real ustar =
        m2 * kappa / (std::log(1.5_rt * dz / z0) - non_neutral_neighbour);
    Dxz +=
        -ustar * ustar * ux1r /
        (amr_wind::constants::EPS + std::sqrt((ux1r * ux1r) + (uy1r * uy1r))) /
        dz;
    Dyz +=
        -ustar * ustar * uy1r /
        (amr_wind::constants::EPS + std::sqrt((ux1r * ux1r) + (uy1r * uy1r))) /
        dz;
    return ustar;
}

// Implementation comes from MOSD approach in boundary_conditions/wall_models/
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void form_drag_calculations(
    amrex::Real& Dxz,
    amrex::Real& Dyz,
    const int i,
    const int j,
    const int k,
    amrex::Array4<amrex::Real const> const& phi,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
    const amrex::Real ux1r,
    const amrex::Real uy1r)
{
    // phi = eta - z, so eta derivatives in x and y can be calculate with phi
    amrex::Real n_x, n_y, n_z;
    amr_wind::multiphase::youngs_finite_difference_normal(
        i, j, k, phi, n_x, n_y, n_z);
    // factor of 32 has to do with finite differences, number of points used
    // negative to make normals point away from waves into air
    const amrex::Real dx_eta_wave = -n_x / 32.0_rt / dx[0];
    const amrex::Real dy_eta_wave = -n_y / 32.0_rt / dx[1];
    const amrex::Real grad_eta_wave =
        std::sqrt((dx_eta_wave * dx_eta_wave) + (dy_eta_wave * dy_eta_wave));
    n_x = dx_eta_wave / grad_eta_wave;
    n_y = dy_eta_wave / grad_eta_wave;

    // Relative velocity while considering interface normal
    const amrex::Real ur_mag =
        std::sqrt((ux1r * ux1r * n_x * n_x) + (uy1r * uy1r * n_y * n_y));
    // Heaviside function changes behavior for velocity surplus/deficit
    const amrex::Real Heavi_arg = ((ux1r * dx_eta_wave) + (uy1r * dy_eta_wave));
    const amrex::Real Heavi =
        (Heavi_arg + std::abs(Heavi_arg)) / (2.0_rt * Heavi_arg);

    // Stress in each direction
    const amrex::Real tau_xz =
        (1.0_rt / std::numbers::pi_v<amrex::Real>)*ur_mag * ur_mag *
        grad_eta_wave * grad_eta_wave * Heavi * n_x;
    const amrex::Real tau_yz =
        (1.0_rt / std::numbers::pi_v<amrex::Real>)*ur_mag * ur_mag *
        grad_eta_wave * grad_eta_wave * Heavi * n_y;

    // Drag terms from waves added to log law
    Dxz += -tau_xz / dx[2];
    Dyz += -tau_yz / dx[2];
}
} // namespace

namespace amr_wind::pde::icns {

DragForcing::DragForcing(const CFDSim& sim)
    : m_time(sim.time())
    , m_sim(sim)
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
{
    amrex::ParmParse pp("DragForcing");
    pp.query("drag_coefficient", m_drag_coefficient);
    pp.query("sponge_strength", m_sponge_strength);
    pp.query("sponge_density", m_sponge_density);
    pp.query("sponge_distance_west", m_sponge_distance_west);
    pp.query("sponge_distance_east", m_sponge_distance_east);
    pp.query("sponge_distance_south", m_sponge_distance_south);
    pp.query("sponge_distance_north", m_sponge_distance_north);
    pp.query("sponge_west", m_sponge_west);
    pp.query("sponge_east", m_sponge_east);
    pp.query("sponge_south", m_sponge_south);
    pp.query("sponge_north", m_sponge_north);
    pp.query("is_laminar", m_is_laminar);
    const auto& phy_mgr = m_sim.physics_manager();
    if (phy_mgr.contains("ABL")) {
        const auto& abl = m_sim.physics_manager().get<amr_wind::ABL>();
        const auto& fa_velocity = abl.abl_statistics().vel_profile_coarse();
        m_device_vel_ht.resize(fa_velocity.line_centroids().size());
        m_device_vel_vals.resize(fa_velocity.line_average().size());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, fa_velocity.line_centroids().begin(),
            fa_velocity.line_centroids().end(), m_device_vel_ht.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, fa_velocity.line_average().begin(),
            fa_velocity.line_average().end(), m_device_vel_vals.begin());
    } else {
        m_sponge_strength = 0.0_rt;
    }
    if (phy_mgr.contains("OceanWaves") && !sim.repo().field_exists("vof")) {
        const auto terrain_phys =
            m_sim.physics_manager().get<amr_wind::terraindrag::TerrainDrag>();
        const auto target_vel_name = terrain_phys.wave_velocity_field_name();
        m_target_vel = &sim.repo().get_field(target_vel_name);
        const auto target_levelset_name =
            terrain_phys.wave_negative_elevation_name();
        m_target_levelset = &sim.repo().get_field(target_levelset_name);
        m_terrain_is_waves = true;
        // Inviscid form drag model can help when waves are smaller than cells,
        // i.e., too small to be resolved with cell blanking
        pp.query("wave_model_inviscid_form_drag", m_apply_MOSD);
    }
    amrex::ParmParse pp_abl("ABL");
    pp_abl.query("wall_het_model", m_wall_het_model);
    pp_abl.query("monin_obukhov_length", m_monin_obukhov_length);
    pp_abl.query("kappa", m_kappa);
    pp_abl.query("mo_gamma_m", m_gamma_m);
    pp_abl.query("mo_beta_m", m_beta_m);
}

DragForcing::~DragForcing() = default;

void DragForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& vel =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_array(mfi);
    const bool is_terrain =
        this->m_sim.repo().int_field_exists("terrain_blank");
    if (!is_terrain) {
        amrex::Abort("Need terrain blanking variable to use this source term");
    }
    auto* const m_terrain_blank =
        &this->m_sim.repo().get_int_field("terrain_blank");
    const auto& blank = (*m_terrain_blank)(lev).const_array(mfi);
    auto* const m_terrain_drag =
        &this->m_sim.repo().get_int_field("terrain_drag");
    const auto& drag = (*m_terrain_drag)(lev).const_array(mfi);
    auto* const m_terrainz0 = &this->m_sim.repo().get_field("terrainz0");
    const auto& terrainz0 = (*m_terrainz0)(lev).const_array(mfi);

    const bool is_waves = m_terrain_is_waves;
    const bool model_form_drag = m_apply_MOSD;
    const auto& target_vel_arr = is_waves
                                     ? (*m_target_vel)(lev).const_array(mfi)
                                     : amrex::Array4<amrex::Real>();
    const auto& target_lvs_arr =
        is_waves ? (*m_target_levelset)(lev).const_array(mfi)
                 : amrex::Array4<amrex::Real>();

    const auto& geom = m_mesh.Geom(lev);
    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();
    const amrex::Real drag_coefficient = m_drag_coefficient;
    const amrex::Real sponge_strength = m_sponge_strength;
    const amrex::Real sponge_density = m_sponge_density;
    const amrex::Real start_east = prob_hi[0] - m_sponge_distance_east;
    const amrex::Real start_west = prob_lo[0] - m_sponge_distance_west;
    const amrex::Real start_north = prob_hi[1] - m_sponge_distance_north;
    const amrex::Real start_south = prob_lo[1] - m_sponge_distance_south;
    const int sponge_east = m_sponge_east;
    const int sponge_west = m_sponge_west;
    const int sponge_south = m_sponge_south;
    const int sponge_north = m_sponge_north;
    // Copy Data
    const auto* device_vel_ht = m_device_vel_ht.data();
    const auto* device_vel_vals = m_device_vel_vals.data();
    const unsigned vsize = m_device_vel_ht.size();
    const auto& dt = m_time.delta_t();
    const bool is_laminar = m_is_laminar;
    const amrex::Real scale_factor = (dx[2] < 1.0_rt) ? 1.0_rt : 1.0_rt / dx[2];
    const amrex::Real Cd =
        (is_laminar && dx[2] < 1) ? drag_coefficient : drag_coefficient / dx[2];
    const amrex::Real z0_min = 1.0e-4_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real cd_max = 1000.0_rt;

    const amrex::Real non_neutral_neighbour =
        (m_wall_het_model == "mol")
            ? MOData::calc_psi_m(
                  1.5_rt * dx[2] / m_monin_obukhov_length, m_beta_m, m_gamma_m)
            : 0.0_rt;
    const amrex::Real non_neutral_cell =
        (m_wall_het_model == "mol")
            ? MOData::calc_psi_m(
                  0.5_rt * dx[2] / m_monin_obukhov_length, m_beta_m, m_gamma_m)
            : 0.0_rt;
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        const amrex::Real x = prob_lo[0] + ((i + 0.5_rt) * dx[0]);
        const amrex::Real y = prob_lo[1] + ((j + 0.5_rt) * dx[1]);
        const amrex::Real z = prob_lo[2] + ((k + 0.5_rt) * dx[2]);
        amrex::Real xstart_damping = 0.0_rt;
        amrex::Real ystart_damping = 0.0_rt;
        amrex::Real xend_damping = 0.0_rt;
        amrex::Real yend_damping = 0.0_rt;
        amrex::Real xi_end = (x - start_east) / (prob_hi[0] - start_east);
        amrex::Real xi_start = (start_west - x) / (start_west - prob_lo[0]);
        xi_start = sponge_west * amrex::max<amrex::Real>(xi_start, 0.0_rt);
        xi_end = sponge_east * amrex::max<amrex::Real>(xi_end, 0.0_rt);
        xstart_damping = sponge_west * sponge_strength * xi_start * xi_start;
        xend_damping = sponge_east * sponge_strength * xi_end * xi_end;
        amrex::Real yi_end = (y - start_north) / (prob_hi[1] - start_north);
        amrex::Real yi_start = (start_south - y) / (start_south - prob_lo[1]);
        yi_start = sponge_south * amrex::max<amrex::Real>(yi_start, 0.0_rt);
        yi_end = sponge_north * amrex::max<amrex::Real>(yi_end, 0.0_rt);
        ystart_damping = sponge_strength * yi_start * yi_start;
        yend_damping = sponge_strength * yi_end * yi_end;
        const amrex::Real ux1 = vel(i, j, k, 0);
        const amrex::Real uy1 = vel(i, j, k, 1);
        const amrex::Real uz1 = vel(i, j, k, 2);
        const auto idx =
            interp::bisection_search(device_vel_ht, device_vel_ht + vsize, z);
        const amrex::Real spongeVelX =
            (vsize > 0) ? interp::linear_impl(
                              device_vel_ht, device_vel_vals, z, idx, 3, 0)
                        : ux1;
        const amrex::Real spongeVelY =
            (vsize > 0) ? interp::linear_impl(
                              device_vel_ht, device_vel_vals, z, idx, 3, 1)
                        : uy1;
        const amrex::Real spongeVelZ =
            (vsize > 0) ? interp::linear_impl(
                              device_vel_ht, device_vel_vals, z, idx, 3, 2)
                        : uz1;
        amrex::Real Dxz = 0.0_rt;
        amrex::Real Dyz = 0.0_rt;
        amrex::Real bc_forcing_x = 0.0_rt;
        amrex::Real bc_forcing_y = 0.0_rt;
        const amrex::Real m =
            std::sqrt((ux1 * ux1) + (uy1 * uy1) + (uz1 * uz1));
        if (drag(i, j, k) == 1 && (!is_laminar)) {
            // Check if close enough to interface to use current cell or below
            int k_off = -1;
            if (is_waves) {
                const amrex::Real cell_length_2D =
                    std::sqrt((dx[0] * dx[0]) + (dx[2] * dx[2]));
                if (target_lvs_arr(i, j, k) + cell_length_2D >= 0) {
                    // Current cell will be used for wave velocity
                    k_off = 0;
                }
                // Cell below will be used if not (default of -1)
            }
            // Establish wall velocity
            // - estimate wave velocity using target velocity in cells below
            const amrex::Real wall_u =
                !is_waves ? 0.0_rt : target_vel_arr(i, j, k + k_off, 0);
            const amrex::Real wall_v =
                !is_waves ? 0.0_rt : target_vel_arr(i, j, k + k_off, 1);
            // Relative velocities for calculating shear
            const amrex::Real ux1r = ux1 - wall_u;
            const amrex::Real uy1r = uy1 - wall_v;
            const amrex::Real ux2r = vel(i, j, k + 1, 0) - wall_u;
            const amrex::Real uy2r = vel(i, j, k + 1, 1) - wall_v;
            const amrex::Real z0 =
                amrex::max<amrex::Real>(terrainz0(i, j, k), z0_min);
            const amrex::Real ustar = viscous_drag_calculations(
                Dxz, Dyz, ux1r, uy1r, ux2r, uy2r, z0, dx[2], kappa,
                non_neutral_neighbour);
            if (model_form_drag) {
                form_drag_calculations(
                    Dxz, Dyz, i, j, k, target_lvs_arr, dx, ux1r, uy1r);
            }
            const amrex::Real uTarget =
                ustar / kappa *
                (std::log(0.5_rt * dx[2] / z0) - non_neutral_cell);
            const amrex::Real uxTarget =
                uTarget * ux2r /
                (amr_wind::constants::EPS +
                 std::sqrt((ux2r * ux2r) + (uy2r * uy2r)));
            const amrex::Real uyTarget =
                uTarget * uy2r /
                (amr_wind::constants::EPS +
                 std::sqrt((ux2r * ux2r) + (uy2r * uy2r)));
            // BC forcing pushes nonrelative velocity toward target velocity
            bc_forcing_x = -(uxTarget - ux1) / dt;
            bc_forcing_y = -(uyTarget - uy1) / dt;
        }
        // Target velocity intended for within terrain
        amrex::Real target_u = 0.0_rt;
        amrex::Real target_v = 0.0_rt;
        amrex::Real target_w = 0.0_rt;
        if (is_waves) {
            target_u = target_vel_arr(i, j, k, 0);
            target_v = target_vel_arr(i, j, k, 1);
            target_w = target_vel_arr(i, j, k, 2);
        }

        const amrex::Real CdM = amrex::min<amrex::Real>(
            Cd / (m + amr_wind::constants::EPS), cd_max / scale_factor);
        src_term(i, j, k, 0) -=
            ((CdM * m * (ux1 - target_u) * blank(i, j, k)) +
             (Dxz * drag(i, j, k)) + (bc_forcing_x * drag(i, j, k)) +
             ((xstart_damping + xend_damping + ystart_damping + yend_damping) *
              (ux1 - sponge_density * spongeVelX)));
        src_term(i, j, k, 1) -=
            ((CdM * m * (uy1 - target_v) * blank(i, j, k)) +
             (Dyz * drag(i, j, k)) + (bc_forcing_y * drag(i, j, k)) +
             ((xstart_damping + xend_damping + ystart_damping + yend_damping) *
              (uy1 - sponge_density * spongeVelY)));
        src_term(i, j, k, 2) -=
            ((CdM * m * (uz1 - target_w) * blank(i, j, k)) +
             ((xstart_damping + xend_damping + ystart_damping + yend_damping) *
              (uz1 - sponge_density * spongeVelZ)));
    });
}

} // namespace amr_wind::pde::icns
