#include "AMReX_Orientation.H"
#include "AMReX_Gpu.H"
#include "src/equation_systems/tke/source_terms/KransAxell.H"
#include "src/CFDSim.H"
#include "src/turbulence/TurbulenceModel.H"
#include "src/wind_energy/MOData.H"
#include "src/utilities/linear_interpolation.H"
#include "src/utilities/constants.H"
#include "src/utilities/math_ops.H"

using namespace amrex::literals;
namespace kynema_sgf::pde::tke {

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
    const std::string model_name = sim.turbulence_model().model_name();
    AMREX_ALWAYS_ASSERT(
        model_name == "KLAxell" || model_name == "KLAxellSeparation");
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

    amrex::ParmParse pp_incflo("incflo");
    pp_incflo.queryarr("gravity", m_gravity);

    amrex::ParmParse pp_drag("DragForcing");
    pp_drag.query("bc_forcing_time_factor", m_forcing_time_factor);
    pp_drag.query("sponge_west", m_sponge_west);
    pp_drag.query("sponge_east", m_sponge_east);
    pp_drag.query("sponge_south", m_sponge_south);
    pp_drag.query("sponge_north", m_sponge_north);
    if (m_sponge_west) {
        pp_drag.get("sponge_distance_west", m_sponge_distance_west);
    }
    if (m_sponge_east) {
        pp_drag.get("sponge_distance_east", m_sponge_distance_east);
    }
    if (m_sponge_south) {
        pp_drag.get("sponge_distance_south", m_sponge_distance_south);
    }
    if (m_sponge_north) {
        pp_drag.get("sponge_distance_north", m_sponge_distance_north);
    }

    if (model_name == "KLAxellSeparation") {
        amrex::ParmParse pp_separation("KLAxellSeparation");
        pp_separation.query("production_cap", m_production_cap);
        pp_separation.query("destruction_boost", m_destruction_boost);
        pp_separation.query("implicit_dissipation", m_implicit_dissipation);
        if (m_implicit_dissipation) {
            // The dissipation moves to the diagonal of the TKE diffusion
            // solve, which the explicit diffusion type does not have. Same
            // default diffusion type (Crank-Nicolson) as in init.cpp.
            int diffusion_type = 1;
            pp_incflo.query("diffusion_type", diffusion_type);
            if (diffusion_type == 2) {
                m_implicit_factor = 1.0_rt;
            } else if (diffusion_type == 1) {
                m_implicit_factor = 0.5_rt;
            } else {
                amrex::Abort(
                    "KLAxellSeparation.implicit_dissipation requires "
                    "incflo.diffusion_type = 1 (Crank-Nicolson) or 2 "
                    "(implicit)");
            }
        }
        // The model declares the gate field when the sensor is on and its
        // relaxation time is positive (the default)
        m_relaxed_gate = sim.repo().field_exists("separation_gate");
    }
}

KransAxell::~KransAxell() = default;

void KransAxell::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    if (!m_ref_theta_scratch) {
        m_ref_theta_scratch = m_sim.repo().create_scratch_field(1, 0);
    }
    m_transport.ref_theta_fill(lev, (*m_ref_theta_scratch)(lev));
    amrex::Gpu::streamSynchronize();

    const auto& geom = m_mesh.Geom(lev);
    const auto& problo = geom.ProbLoArray();
    const auto& probhi = geom.ProbHiArray();
    const auto& dx = geom.CellSizeArray();
    const auto& dt = m_time.delta_t();
    const amrex::Real heat_flux = m_heat_flux;
    const amrex::Real Cmu = m_Cmu;
    const amrex::Real kappa = m_kappa;
    const amrex::Real time_factor = m_forcing_time_factor;
    const amrex::Real z0 = m_z0;
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

    auto const& src_arrs = src_term.arrays();
    auto const& vel_arrs =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_arrays();
    auto const& tlscale_arrs = m_turb_lscale(lev).const_arrays();
    auto const& shear_prod_arrs = m_shear_prod(lev).const_arrays();
    auto const& buoy_prod_arrs = m_buoy_prod(lev).const_arrays();
    auto const& dissip_arrs = m_dissip(lev).arrays();
    auto const& tke_arrs = m_tke(lev).const_arrays();
    auto const& ref_theta_arrs = (*m_ref_theta_scratch)(lev).const_arrays();

    // Main kernel
    amrex::ParallelFor(
        src_term, amrex::IntVect(0), 1,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
            const auto& vel = vel_arrs[nbx];
            const auto& tlscale_arr = tlscale_arrs[nbx];
            const auto& shear_prod_arr = shear_prod_arrs[nbx];
            const auto& buoy_prod_arr = buoy_prod_arrs[nbx];
            const auto& dissip_arr = dissip_arrs[nbx];
            const auto& tke_arr = tke_arrs[nbx];
            const auto& ref_theta_arr = ref_theta_arrs[nbx];

            amrex::Real bcforcing = 0;
            const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
            if (k == 0) {
                const amrex::Real ux = vel(i, j, k + 1, 0);
                const amrex::Real uy = vel(i, j, k + 1, 1);
                const amrex::Real m = std::sqrt((ux * ux) + (uy * uy));
                const amrex::Real ustar =
                    m * kappa / (std::log(3.0_rt * z / z0) - psi_m);
                const amrex::Real T0 = ref_theta_arr(i, j, k);
                const amrex::Real hf = std::abs(gravity[2]) / T0 * heat_flux;
                const amrex::Real rans_b = amrex::max<amrex::Real>(hf, 0.0_rt) *
                                           kappa * z / utils::powi(Cmu, 3);
                const amrex::Real tke_exact = std::pow(
                    ustar * ustar * ustar / (Cmu * Cmu * Cmu) + rans_b,
                    2.0_rt / 3.0_rt);
                bcforcing = (tke_exact - tke_arr(i, j, k)) / (time_factor * dt);
            }
            amrex::Real ref_tke = tke_arr(i, j, k);
            if (z > sponge_start) {
                ref_tke = (vsize > 0)
                              ? interp::linear(
                                    wind_heights_d, wind_heights_d + vsize,
                                    tke_values_d, z)
                              : tke_arr(i, j, k, 0);
            }
            const amrex::Real sponge_forcing =
                1.0_rt / dt * (tke_arr(i, j, k) - ref_tke);
            dissip_arr(i, j, k) =
                utils::powi(Cmu, 3) * std::pow(tke_arr(i, j, k), 1.5_rt) /
                (tlscale_arr(i, j, k) + kynema_sgf::constants::EPS);
            src_arrs[nbx](i, j, k) +=
                shear_prod_arr(i, j, k) + buoy_prod_arr(i, j, k) -
                dissip_arr(i, j, k) -
                ((1.0_rt - static_cast<int>(has_terrain)) *
                 (sponge_forcing - bcforcing));
        });

    if (has_terrain) {
        const amrex::Real z0_min = 1.0e-4_rt;
        const auto& terrain_blank =
            this->m_sim.repo().get_int_field("terrain_blank");
        const auto& terrain_drag =
            this->m_sim.repo().get_int_field("terrain_drag");
        auto& terrain_height = this->m_sim.repo().get_field("terrain_height");
        auto& terrainz0 = this->m_sim.repo().get_field("terrainz0");

        auto const& blank_arrs = terrain_blank(lev).const_arrays();
        auto const& drag_arrs = terrain_drag(lev).const_arrays();
        auto const& terrain_height_arrs = terrain_height(lev).const_arrays();
        auto const& terrainz0_arrs = terrainz0(lev).const_arrays();

        amrex::ParallelFor(
            src_term, amrex::IntVect(0), 1,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
                const auto& vel = vel_arrs[nbx];
                const auto& tke_arr = tke_arrs[nbx];
                const auto& blank_arr = blank_arrs[nbx];
                const auto& drag_arr = drag_arrs[nbx];
                const auto& terrain_height_arr = terrain_height_arrs[nbx];
                const auto& terrainz0_arr = terrainz0_arrs[nbx];
                const auto& ref_theta_arr = ref_theta_arrs[nbx];

                const amrex::Real cell_z0 =
                    (drag_arr(i, j, k) *
                     amrex::max<amrex::Real>(terrainz0_arr(i, j, k), z0_min)) +
                    ((1.0_rt - drag_arr(i, j, k)) * z0);
                amrex::Real ux = vel(i, j, k + 1, 0);
                amrex::Real uy = vel(i, j, k + 1, 1);
                amrex::Real z = 0.5_rt * dx[2];
                amrex::Real m = std::sqrt((ux * ux) + (uy * uy));
                const amrex::Real ustar =
                    m * kappa / (std::log(3.0_rt * z / cell_z0) - psi_m);
                const amrex::Real T0 = ref_theta_arr(i, j, k);
                const amrex::Real hf = std::abs(gravity[2]) / T0 * heat_flux;
                const amrex::Real rans_b = amrex::max<amrex::Real>(hf, 0.0_rt) *
                                           kappa * z / utils::powi(Cmu, 3);
                const amrex::Real tke_exact = std::pow(
                    ustar * ustar * ustar / (Cmu * Cmu * Cmu) + rans_b,
                    2.0_rt / 3.0_rt);
                const amrex::Real terrainforcing =
                    (tke_exact - tke_arr(i, j, k)) / (time_factor * dt);
                amrex::Real bcforcing = 0.0_rt;
                if (k == 0) {
                    bcforcing = (1 - blank_arr(i, j, k)) * terrainforcing;
                }
                ux = vel(i, j, k, 0);
                uy = vel(i, j, k, 1);
                const amrex::Real uz = vel(i, j, k, 2);
                m = std::sqrt((ux * ux) + (uy * uy) + (uz * uz));
                const amrex::Real Cd = amrex::min<amrex::Real>(
                    10.0_rt / (dx[2] * m + kynema_sgf::constants::EPS),
                    100.0_rt / dx[2]);
                const amrex::Real dragforcing = -Cd * m * tke_arr(i, j, k, 0);
                z = amrex::max<amrex::Real>(
                    problo[2] + ((k + 0.5_rt) * dx[2]) -
                        terrain_height_arr(i, j, k),
                    0.5_rt * dx[2]);
                amrex::Real ref_tke = tke_arr(i, j, k);
                if (z > sponge_start) {
                    ref_tke = (vsize > 0)
                                  ? interp::linear(
                                        wind_heights_d, wind_heights_d + vsize,
                                        tke_values_d, z)
                                  : tke_arr(i, j, k, 0);
                }
                const amrex::Real sponge_forcing =
                    1.0_rt / dt * (tke_arr(i, j, k) - ref_tke);
                src_arrs[nbx](i, j, k) =
                    ((1.0_rt - blank_arr(i, j, k)) * src_arrs[nbx](i, j, k)) +
                    (drag_arr(i, j, k) * terrainforcing) +
                    (blank_arr(i, j, k) * dragforcing) -
                    (static_cast<int>(has_terrain) *
                     (sponge_forcing - bcforcing));
            });

        if (m_horizontal_sponge) {
            const amrex::Real sponge_strength = m_sponge_strength;
            const amrex::Real start_east = probhi[0] - m_sponge_distance_east;
            const amrex::Real start_west = problo[0] - m_sponge_distance_west;
            const amrex::Real start_north = probhi[1] - m_sponge_distance_north;
            const amrex::Real start_south = problo[1] - m_sponge_distance_south;
            const amrex::Real sdist_east = m_sponge_distance_east;
            const amrex::Real sdist_west = m_sponge_distance_west;
            const amrex::Real sdist_north = m_sponge_distance_north;
            const amrex::Real sdist_south = m_sponge_distance_south;
            const auto sponge_east = static_cast<int>(m_sponge_east);
            const auto sponge_west = static_cast<int>(m_sponge_west);
            const auto sponge_south = static_cast<int>(m_sponge_south);
            const auto sponge_north = static_cast<int>(m_sponge_north);

            amrex::ParallelFor(
                src_term, amrex::IntVect(0), 1,
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
                    const auto& tke_arr = tke_arrs[nbx];
                    const amrex::Real x = problo[0] + ((i + 0.5_rt) * dx[0]);
                    const amrex::Real y = problo[1] + ((j + 0.5_rt) * dx[1]);
                    const amrex::Real z = problo[2] + ((k + 0.5_rt) * dx[2]);
                    amrex::Real xi_end =
                        (std::abs(sdist_east) > kynema_sgf::constants::EPS)
                            ? (x - start_east) / (sdist_east)
                            : 0.0_rt;
                    amrex::Real xi_start =
                        (std::abs(sdist_west) > kynema_sgf::constants::EPS)
                            ? (start_west - x) / (-sdist_west)
                            : 0.0_rt;
                    xi_start =
                        sponge_west * amrex::max<amrex::Real>(xi_start, 0.0_rt);
                    xi_end =
                        sponge_east * amrex::max<amrex::Real>(xi_end, 0.0_rt);
                    xi_start /= (xi_start + kynema_sgf::constants::EPS);
                    xi_end /= (xi_end + kynema_sgf::constants::EPS);
                    const amrex::Real xstart_damping =
                        sponge_strength * xi_start * xi_start;
                    const amrex::Real xend_damping =
                        sponge_strength * xi_end * xi_end;
                    amrex::Real yi_end =
                        (std::abs(sdist_north) > kynema_sgf::constants::EPS)
                            ? (y - start_north) / (sdist_north)
                            : 0.0_rt;
                    amrex::Real yi_start =
                        (std::abs(sdist_south) > kynema_sgf::constants::EPS)
                            ? (start_south - y) / (-sdist_south)
                            : 0.0_rt;
                    yi_start = sponge_south *
                               amrex::max<amrex::Real>(yi_start, 0.0_rt);
                    yi_end =
                        sponge_north * amrex::max<amrex::Real>(yi_end, 0.0_rt);
                    yi_start /= (yi_start + kynema_sgf::constants::EPS);
                    yi_end /= (yi_end + kynema_sgf::constants::EPS);
                    const amrex::Real ystart_damping =
                        sponge_strength * yi_start * yi_start;
                    const amrex::Real yend_damping =
                        sponge_strength * yi_end * yi_end;
                    const amrex::Real ref_tke =
                        (vsize > 0)
                            ? interp::linear(
                                  wind_heights_d, wind_heights_d + vsize,
                                  tke_values_d, z)
                            : tke_arr(i, j, k, 0);
                    const amrex::Real damping_sum =
                        (xstart_damping + xend_damping + ystart_damping +
                         yend_damping + kynema_sgf::constants::EPS);
                    const amrex::Real sponge_forcing =
                        (xstart_damping + xend_damping + ystart_damping +
                         yend_damping) /
                        (damping_sum * dt) * (tke_arr(i, j, k) - ref_tke);
                    src_arrs[nbx](i, j, k, 0) -= sponge_forcing;
                });
        }
    }
    if (m_production_cap && m_relaxed_gate) {
        production_cap_relaxed(lev, src_term);
    } else if (m_production_cap) {
        production_cap(lev, src_term);
    }
    if (m_destruction_boost && m_relaxed_gate) {
        destruction_boost_relaxed(lev, src_term);
    } else if (m_destruction_boost) {
        destruction_boost(lev, src_term);
    }
    if (m_implicit_dissipation && (fstate == FieldState::Old)) {
        // The old-time source is the Godunov forcing of the face states. It
        // keeps the full dissipation, so the advective fluxes are those of the
        // explicit treatment and the steady state stays unchanged; the
        // diagonal is cleared so that a solve with this source stays explicit.
        m_sim.repo().get_field("tke_lhs_src_term")(lev).setVal(0.0_rt);
    } else if (m_implicit_dissipation) {
        implicit_dissipation(lev, fstate, src_term);
    }
    amrex::Gpu::streamSynchronize();
}

void KransAxell::production_cap(const int lev, amrex::MultiFab& src_term) const
{
    const auto coeffs = m_sim.turbulence_model().model_coeffs();
    const amrex::Real ratio = coeffs.at("production_cap_ratio");
    const amrex::Real threshold = coeffs.at("sensor_threshold");

    auto const& src_arrs = src_term.arrays();
    auto const& shear_prod_arrs = m_shear_prod(lev).const_arrays();
    auto const& dissip_arrs = m_dissip(lev).const_arrays();
    auto const& sensor_arrs =
        m_sim.repo().get_field("pressure_gradient_sensor")(lev).const_arrays();

    // Remove the shear production above ratio * dissipation where the sensor
    // fires, with the ramp gate of the realizable Cmu limiter. The gate is
    // exactly 0 below the threshold, so those cells keep the KLAxell source.
    amrex::ParallelFor(
        src_term, amrex::IntVect(0), 1,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
            const amrex::Real gate = amrex::min<amrex::Real>(
                amrex::max<amrex::Real>(
                    (sensor_arrs[nbx](i, j, k) - threshold) / threshold,
                    0.0_rt),
                1.0_rt);
            src_arrs[nbx](i, j, k) -=
                gate * amrex::max<amrex::Real>(
                           shear_prod_arrs[nbx](i, j, k) -
                               (ratio * dissip_arrs[nbx](i, j, k)),
                           0.0_rt);
        });
}

void KransAxell::production_cap_relaxed(
    const int lev, amrex::MultiFab& src_term) const
{
    const amrex::Real ratio =
        m_sim.turbulence_model().model_coeffs().at("production_cap_ratio");

    auto const& src_arrs = src_term.arrays();
    auto const& shear_prod_arrs = m_shear_prod(lev).const_arrays();
    auto const& dissip_arrs = m_dissip(lev).const_arrays();
    auto const& gate_arrs =
        m_sim.repo().get_field("separation_gate")(lev).const_arrays();

    // Same cap as production_cap with the stored, time-relaxed gate
    amrex::ParallelFor(
        src_term, amrex::IntVect(0), 1,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
            src_arrs[nbx](i, j, k) -=
                gate_arrs[nbx](i, j, k) *
                amrex::max<amrex::Real>(
                    shear_prod_arrs[nbx](i, j, k) -
                        (ratio * dissip_arrs[nbx](i, j, k)),
                    0.0_rt);
        });
}

void KransAxell::destruction_boost(
    const int lev, amrex::MultiFab& src_term) const
{
    const auto coeffs = m_sim.turbulence_model().model_coeffs();
    const amrex::Real extra = coeffs.at("destruction_boost_factor") - 1.0_rt;
    const amrex::Real threshold = coeffs.at("sensor_threshold");

    auto const& src_arrs = src_term.arrays();
    auto const& dissip_arrs = m_dissip(lev).const_arrays();
    auto const& sensor_arrs =
        m_sim.repo().get_field("pressure_gradient_sensor")(lev).const_arrays();

    // Remove (c_d - 1) eps where the sensor fires, with the ramp gate of the
    // realizable Cmu limiter. The gate is exactly 0 below the threshold, so
    // those cells keep the KLAxell source.
    amrex::ParallelFor(
        src_term, amrex::IntVect(0), 1,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
            const amrex::Real gate = amrex::min<amrex::Real>(
                amrex::max<amrex::Real>(
                    (sensor_arrs[nbx](i, j, k) - threshold) / threshold,
                    0.0_rt),
                1.0_rt);
            src_arrs[nbx](i, j, k) -= gate * extra * dissip_arrs[nbx](i, j, k);
        });
}

void KransAxell::destruction_boost_relaxed(
    const int lev, amrex::MultiFab& src_term) const
{
    const amrex::Real extra =
        m_sim.turbulence_model().model_coeffs().at("destruction_boost_factor") -
        1.0_rt;

    auto const& src_arrs = src_term.arrays();
    auto const& dissip_arrs = m_dissip(lev).const_arrays();
    auto const& gate_arrs =
        m_sim.repo().get_field("separation_gate")(lev).const_arrays();

    // Same boost as destruction_boost with the stored, time-relaxed gate
    amrex::ParallelFor(
        src_term, amrex::IntVect(0), 1,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
            src_arrs[nbx](i, j, k) -=
                gate_arrs[nbx](i, j, k) * extra * dissip_arrs[nbx](i, j, k);
        });
}

void KransAxell::implicit_dissipation(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    const amrex::Real dt = m_time.delta_t();
    const amrex::Real Cmu3 = utils::powi(m_Cmu, 3);
    const amrex::Real factor = m_implicit_factor;

    auto const& src_arrs = src_term.arrays();
    auto const& dissip_arrs = m_dissip(lev).const_arrays();
    auto const& tke_arrs = m_tke(lev).const_arrays();
    auto const& tlscale_arrs = m_turb_lscale(lev).const_arrays();
    auto const& lhs_arrs =
        m_sim.repo().get_field("tke_lhs_src_term")(lev).arrays();
    // Same density state as the multiplication of the source in SrcTermOp
    auto const& rho_arrs = m_sim.repo()
                               .get_field("density")
                               .state(field_impl::phi_state(fstate))(lev)
                               .const_arrays();

    // With eps = (Cmu^3 sqrt(k) / L) k, factor * eps leaves the explicit
    // source and factor * rho * dt * Cmu^3 sqrt(k) / L enters the diagonal of
    // the diffusion solve, with the length scale floored as in the dissipation
    // of the main kernel. For k^{n+1} = k^n the diagonal term equals the
    // removed source, so steady states are unchanged.
    if (m_sim.repo().int_field_exists("terrain_blank")) {
        // Blanked terrain cells have no dissipation in the source
        auto const& blank_arrs =
            m_sim.repo().get_int_field("terrain_blank")(lev).const_arrays();
        amrex::ParallelFor(
            src_term, amrex::IntVect(0), 1,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
                const amrex::Real fluid =
                    1.0_rt - static_cast<amrex::Real>(blank_arrs[nbx](i, j, k));
                src_arrs[nbx](i, j, k) +=
                    fluid * factor * dissip_arrs[nbx](i, j, k);
                lhs_arrs[nbx](i, j, k) =
                    fluid * factor * rho_arrs[nbx](i, j, k) * dt * Cmu3 *
                    std::sqrt(
                        amrex::max<amrex::Real>(
                            tke_arrs[nbx](i, j, k), 0.0_rt)) /
                    (tlscale_arrs[nbx](i, j, k) + kynema_sgf::constants::EPS);
            });
    } else {
        amrex::ParallelFor(
            src_term, amrex::IntVect(0), 1,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int) {
                src_arrs[nbx](i, j, k) += factor * dissip_arrs[nbx](i, j, k);
                lhs_arrs[nbx](i, j, k) =
                    factor * rho_arrs[nbx](i, j, k) * dt * Cmu3 *
                    std::sqrt(
                        amrex::max<amrex::Real>(
                            tke_arrs[nbx](i, j, k), 0.0_rt)) /
                    (tlscale_arrs[nbx](i, j, k) + kynema_sgf::constants::EPS);
            });
    }
}

} // namespace kynema_sgf::pde::tke
