#include <AMReX_Orientation.H>

#include "amr-wind/equation_systems/tke/source_terms/KransAxell.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/turbulence/TurbulenceModel.H"
#include "amr-wind/wind_energy/MOData.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "amr-wind/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;
namespace amr_wind::pde::tke {

KransAxell::KransAxell(const CFDSim& sim)
    : m_turb_lscale(sim.repo().get_field("turb_lscale"))
    , m_shear_prod(sim.repo().get_field("shear_prod"))
    , m_buoy_prod(sim.repo().get_field("buoy_prod"))
    , m_dissip(sim.repo().get_field("dissipation"))
    , m_tke(sim.repo().get_field("tke"))
    , m_time(sim.time())
    , m_sim(sim)
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_transport(sim.transport_model())
{
    AMREX_ALWAYS_ASSERT(sim.turbulence_model().model_name() == "KLAxell");
    auto coeffs = sim.turbulence_model().model_coeffs();
    amrex::ParmParse pp("ABL");
    pp.query("Cmu", m_Cmu);
    pp.query("kappa", m_kappa);
    pp.query("surface_roughness_z0", m_z0);
    pp.query("surface_temp_flux", m_heat_flux);
    pp.query("meso_sponge_start", m_meso_start);
    pp.query("rans_1dprofile_file", m_1d_rans);
    pp.query("horizontal_sponge_tke", m_horizontal_sponge);
    if (!m_1d_rans.empty()) {
        std::ifstream ransfile(m_1d_rans, std::ios::in);
        if (!ransfile.good()) {
            amrex::Abort("Cannot find 1-D RANS profile file " + m_1d_rans);
        }
        amrex::Real value1, value2, value3, value4, value5;
        while (ransfile >> value1 >> value2 >> value3 >> value4 >> value5) {
            m_wind_heights.push_back(value1);
            m_tke_values.push_back(value5);
        }
        int num_wind_values = static_cast<int>(m_wind_heights.size());
        m_wind_heights_d.resize(num_wind_values);
        m_tke_values_d.resize(num_wind_values);
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_wind_heights.begin(),
            m_wind_heights.end(), m_wind_heights_d.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, m_tke_values.begin(), m_tke_values.end(),
            m_tke_values_d.begin());
    } else {
        amrex::Abort("Cannot find 1-D RANS profile file " + m_1d_rans);
    }
    pp.query("wall_het_model", m_wall_het_model);
    pp.query("monin_obukhov_length", m_monin_obukhov_length);
    pp.query("mo_gamma_m", m_gamma_m);
    pp.query("mo_beta_m", m_beta_m);
    {
        amrex::ParmParse pp_incflow("incflo");
        pp_incflow.queryarr("gravity", m_gravity);
    }

    amrex::ParmParse pp_drag("DragForcing");
    pp_drag.query("sponge_strength", m_sponge_strength);
    pp_drag.query("sponge_distance_west", m_sponge_distance_west);
    pp_drag.query("sponge_distance_east", m_sponge_distance_east);
    pp_drag.query("sponge_distance_south", m_sponge_distance_south);
    pp_drag.query("sponge_distance_north", m_sponge_distance_north);
    pp_drag.query("sponge_west", m_sponge_west);
    pp_drag.query("sponge_east", m_sponge_east);
    pp_drag.query("sponge_south", m_sponge_south);
    pp_drag.query("sponge_north", m_sponge_north);
}

KransAxell::~KransAxell() = default;

void KransAxell::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState fstate,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& vel =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_array(mfi);
    const auto& tlscale_arr = (this->m_turb_lscale)(lev).array(mfi);
    const auto& shear_prod_arr = (this->m_shear_prod)(lev).array(mfi);
    const auto& buoy_prod_arr = (this->m_buoy_prod)(lev).array(mfi);
    const auto& dissip_arr = (this->m_dissip)(lev).array(mfi);
    const auto& tke_arr = m_tke(lev).array(mfi);
    amrex::FArrayBox ref_theta_fab(bx, 1, amrex::The_Async_Arena());
    amrex::Array4<amrex::Real> const& ref_theta_arr = ref_theta_fab.array();
    m_transport.ref_theta_impl(lev, mfi, bx, ref_theta_arr);
    const auto& geom = m_mesh.Geom(lev);
    const auto& problo = m_mesh.Geom(lev).ProbLoArray();
    const auto& probhi = m_mesh.Geom(lev).ProbHiArray();
    const auto& dx = geom.CellSizeArray();
    const auto& dt = m_time.delta_t();
    const amrex::Real heat_flux = m_heat_flux;
    const amrex::Real Cmu = m_Cmu;
    const amrex::Real kappa = m_kappa;
    amrex::Real z0 = m_z0;
    const bool has_terrain =
        this->m_sim.repo().int_field_exists("terrain_blank");
    const amrex::Real sponge_start = m_meso_start;
    const auto vsize = m_wind_heights_d.size();
    const auto* wind_heights_d = m_wind_heights_d.data();
    const auto* tke_values_d = m_tke_values_d.data();
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};
    amrex::Real psi_m = 0.0_rt;
    if (m_wall_het_model == "mol") {
        psi_m = MOData::calc_psi_m(
            1.5_rt * dx[2] / m_monin_obukhov_length, m_beta_m, m_gamma_m);
    }
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        amrex::Real bcforcing = 0;
        const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
        if (k == 0) {
            const amrex::Real ux = vel(i, j, k + 1, 0);
            const amrex::Real uy = vel(i, j, k + 1, 1);
            const amrex::Real m = std::sqrt((ux * ux) + (uy * uy));
            const amrex::Real ustar =
                m * kappa / (std::log(3 * z / z0) - psi_m);
            const amrex::Real T0 = ref_theta_arr(i, j, k);
            const amrex::Real hf = std::abs(gravity[2]) / T0 * heat_flux;
            const amrex::Real rans_b = std::pow(
                amrex::max<amrex::Real>(hf, 0.0_rt) * kappa * z /
                    std::pow(Cmu, 3.0_rt),
                (2.0_rt / 3.0_rt));
            bcforcing =
                (ustar * ustar / (Cmu * Cmu) + rans_b - tke_arr(i, j, k)) / dt;
        }
        amrex::Real ref_tke = tke_arr(i, j, k);
        const amrex::Real zi = amrex::max<amrex::Real>(
            (z - sponge_start) / (probhi[2] - sponge_start), 0.0_rt);
        if (zi > 0) {
            ref_tke = (vsize > 0) ? interp::linear(
                                        wind_heights_d, wind_heights_d + vsize,
                                        tke_values_d, z)
                                  : tke_arr(i, j, k, 0);
        }
        const amrex::Real sponge_forcing =
            1.0_rt / dt * (tke_arr(i, j, k) - ref_tke);
        dissip_arr(i, j, k) = std::pow(Cmu, 3.0_rt) *
                              std::pow(tke_arr(i, j, k), 1.5_rt) /
                              (tlscale_arr(i, j, k) + amr_wind::constants::EPS);
        src_term(i, j, k) += shear_prod_arr(i, j, k) + buoy_prod_arr(i, j, k) -
                             dissip_arr(i, j, k) -
                             ((1.0_rt - static_cast<int>(has_terrain)) *
                              (sponge_forcing - bcforcing));
    });
    if (has_terrain) {
        const amrex::Real z0_min = 1.0e-4_rt;
        const auto* const m_terrain_blank =
            &this->m_sim.repo().get_int_field("terrain_blank");
        const auto* const m_terrain_drag =
            &this->m_sim.repo().get_int_field("terrain_drag");
        auto* const m_terrain_height =
            &this->m_sim.repo().get_field("terrain_height");
        auto* const m_terrainz0 = &this->m_sim.repo().get_field("terrainz0");
        const auto& blank_arr = (*m_terrain_blank)(lev).const_array(mfi);
        const auto& drag_arr = (*m_terrain_drag)(lev).const_array(mfi);
        const auto& terrain_height = (*m_terrain_height)(lev).const_array(mfi);
        const auto& terrainz0 = (*m_terrainz0)(lev).const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            const amrex::Real cell_z0 =
                (drag_arr(i, j, k) *
                 amrex::max<amrex::Real>(terrainz0(i, j, k), z0_min)) +
                ((1.0_rt - drag_arr(i, j, k)) * z0);
            amrex::Real terrainforcing = 0;
            amrex::Real dragforcing = 0;
            amrex::Real ux = vel(i, j, k + 1, 0);
            amrex::Real uy = vel(i, j, k + 1, 1);
            amrex::Real z = 0.5_rt * dx[2];
            amrex::Real m = std::sqrt((ux * ux) + (uy * uy));
            const amrex::Real ustar =
                m * kappa / (std::log(3.0_rt * z / cell_z0) - psi_m);
            const amrex::Real T0 = ref_theta_arr(i, j, k);
            const amrex::Real hf = std::abs(gravity[2]) / T0 * heat_flux;
            const amrex::Real rans_b = std::pow(
                amrex::max<amrex::Real>(hf, 0.0_rt) * kappa * z /
                    std::pow(Cmu, 3.0_rt),
                (2.0_rt / 3.0_rt));
            terrainforcing =
                (ustar * ustar / (Cmu * Cmu) + rans_b - tke_arr(i, j, k)) / dt;
            amrex::Real bcforcing = 0.0_rt;
            if (k == 0) {
                bcforcing = (1 - blank_arr(i, j, k)) * terrainforcing;
            }
            ux = vel(i, j, k, 0);
            uy = vel(i, j, k, 1);
            const amrex::Real uz = vel(i, j, k, 2);
            m = std::sqrt((ux * ux) + (uy * uy) + (uz * uz));
            const amrex::Real Cd = amrex::min<amrex::Real>(
                10.0_rt / (dx[2] * m + amr_wind::constants::EPS),
                100.0_rt / dx[2]);
            dragforcing = -Cd * m * tke_arr(i, j, k, 0);
            z = amrex::max<amrex::Real>(
                problo[2] + ((k + 0.5_rt) * dx[2]) - terrain_height(i, j, k),
                0.5_rt * dx[2]);
            const amrex::Real zi = amrex::max<amrex::Real>(
                (z - sponge_start) / (probhi[2] - sponge_start), 0.0_rt);
            amrex::Real ref_tke = tke_arr(i, j, k);
            if (zi > 0) {
                ref_tke = (vsize > 0)
                              ? interp::linear(
                                    wind_heights_d, wind_heights_d + vsize,
                                    tke_values_d, z)
                              : tke_arr(i, j, k, 0);
            }
            const amrex::Real sponge_forcing =
                1.0_rt / dt * (tke_arr(i, j, k) - ref_tke);
            src_term(i, j, k) =
                ((1.0_rt - blank_arr(i, j, k)) * src_term(i, j, k)) +
                (drag_arr(i, j, k) * terrainforcing) +
                (blank_arr(i, j, k) * dragforcing) -
                (static_cast<int>(has_terrain) * (sponge_forcing - bcforcing));
        });
        if (m_horizontal_sponge) {
            const amrex::Real sponge_strength = m_sponge_strength;
            const amrex::Real start_east = probhi[0] - m_sponge_distance_east;
            const amrex::Real start_west = problo[0] - m_sponge_distance_west;
            const amrex::Real start_north = probhi[1] - m_sponge_distance_north;
            const amrex::Real start_south = problo[1] - m_sponge_distance_south;
            const int sponge_east = m_sponge_east;
            const int sponge_west = m_sponge_west;
            const int sponge_south = m_sponge_south;
            const int sponge_north = m_sponge_north;
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
                const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
                amrex::Real xstart_damping = 0.0_rt;
                amrex::Real ystart_damping = 0.0_rt;
                amrex::Real xend_damping = 0.0_rt;
                amrex::Real yend_damping = 0.0_rt;
                amrex::Real xi_end =
                    (x - start_east) / (probhi[0] - start_east);
                amrex::Real xi_start =
                    (start_west - x) / (start_west - problo[0]);
                xi_start =
                    sponge_west * amrex::max<amrex::Real>(xi_start, 0.0_rt);
                xi_end = sponge_east * amrex::max<amrex::Real>(xi_end, 0.0_rt);
                xi_start /= (xi_start + amr_wind::constants::EPS);
                xi_end /= (xi_end + amr_wind::constants::EPS);
                xstart_damping =
                    sponge_west * sponge_strength * xi_start * xi_start;
                xend_damping = sponge_east * sponge_strength * xi_end * xi_end;
                amrex::Real yi_end =
                    (y - start_north) / (probhi[1] - start_north);
                amrex::Real yi_start =
                    (start_south - y) / (start_south - problo[1]);
                yi_start =
                    sponge_south * amrex::max<amrex::Real>(yi_start, 0.0_rt);
                yi_end = sponge_north * amrex::max<amrex::Real>(yi_end, 0.0_rt);
                yi_start /= (yi_start + amr_wind::constants::EPS);
                yi_end /= (yi_end + amr_wind::constants::EPS);
                ystart_damping = sponge_strength * yi_start * yi_start;
                yend_damping = sponge_strength * yi_end * yi_end;
                const amrex::Real ref_tke =
                    (vsize > 0) ? interp::linear(
                                      wind_heights_d, wind_heights_d + vsize,
                                      tke_values_d, z)
                                : tke_arr(i, j, k, 0);
                const amrex::Real damping_sum =
                    (xstart_damping + xend_damping + ystart_damping +
                     yend_damping + amr_wind::constants::EPS);
                const amrex::Real sponge_forcing =
                    (xstart_damping + xend_damping + ystart_damping +
                     yend_damping) /
                    (damping_sum * dt) * (tke_arr(i, j, k) - ref_tke);
                src_term(i, j, k, 0) -= sponge_forcing;
            });
        }
    }
}

} // namespace amr_wind::pde::tke
