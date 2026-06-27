#include "src/equation_systems/temperature/source_terms/EBDragTempForcing.H"
#include "src/utilities/IOManager.H"
#include "src/wind_energy/MOData.H"
#include "AMReX_ParmParse.H"
#include "AMReX_Gpu.H"
#include "AMReX_Random.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::temperature {

EBDragTempForcing::EBDragTempForcing(const CFDSim& sim)
    : m_time(sim.time())
    , m_sim(sim)
    , m_mesh(sim.mesh())
    , m_velocity(sim.repo().get_field("velocity"))
    , m_temperature(sim.repo().get_field("temperature"))
{
    amrex::ParmParse pp("EBDragTempForcing");
    pp.query("drag_coefficient", m_drag_coefficient);
    pp.query("soil_temperature", m_soil_temperature);
    pp.query("bc_forcing_time_factor", m_forcing_time_factor);
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

    // Parse IDW soil temperature model parameters
    pp.query("use_idw_soil_temp_model", m_use_idw_model);
    if (m_use_idw_model) {
        amrex::Vector<amrex::Real> idw_x_host, idw_y_host, idw_z_host,
            idw_temp_host;
        pp.getarr("idw_x", idw_x_host);
        pp.getarr("idw_y", idw_y_host);
        pp.getarr("idw_z", idw_z_host);
        pp.getarr("idw_temp", idw_temp_host);

        m_idw_num_points = static_cast<int>(idw_x_host.size());

        // Validate that all arrays have the same size
        if ((idw_y_host.size() != m_idw_num_points) ||
            (idw_z_host.size() != m_idw_num_points) ||
            (idw_temp_host.size() != m_idw_num_points)) {
            amrex::Abort(
                "EBDragTempForcing: idw_x, idw_y, idw_z, and idw_temp must "
                "have the same length");
        }

        // Copy to device vectors
        m_idw_x.resize(m_idw_num_points);
        m_idw_y.resize(m_idw_num_points);
        m_idw_z.resize(m_idw_num_points);
        m_idw_temp.resize(m_idw_num_points);

        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, idw_x_host.begin(), idw_x_host.end(),
            m_idw_x.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, idw_y_host.begin(), idw_y_host.end(),
            m_idw_y.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, idw_z_host.begin(), idw_z_host.end(),
            m_idw_z.begin());
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, idw_temp_host.begin(),
            idw_temp_host.end(), m_idw_temp.begin());

        amrex::Print()
            << "EBDragTempForcing: Using IDW soil temperature model with "
            << m_idw_num_points << " points\n";
        // Parse heat flux soil temperature model parameters
        pp.query("use_heatflux_model", m_use_heatflux_model);
        if (m_use_heatflux_model) {
            pp.query("surface_heatflux", m_surface_heatflux);
            pp.query("thermal_conductivity", m_thermal_conductivity);
            pp.query("boundary_layer_height", m_boundary_layer_height);

            amrex::Print()
                << "EBDragTempForcing: Using heat flux soil temperature model\n"
                << "  Surface heat flux: " << m_surface_heatflux << " W/m^2\n"
                << "  Thermal conductivity: " << m_thermal_conductivity
                << " W/(m·K)\n"
                << "  Boundary layer height: " << m_boundary_layer_height
                << " m\n";
        }
    }
}

EBDragTempForcing::~EBDragTempForcing() = default;

void EBDragTempForcing::operator()(
    const int lev, const FieldState fstate, amrex::MultiFab& src_term) const
{
    auto const& src_arrs = src_term.arrays();
    auto const& vel =
        m_velocity.state(field_impl::dof_state(fstate))(lev).const_arrays();
    auto const& temperature =
        m_temperature.state(field_impl::dof_state(fstate))(lev).const_arrays();
    const int is_eb = this->m_sim.repo().field_exists("eb_blank") ? 1 : 0;
    if (is_eb == 0) {
        amrex::Abort("Need EB blanking variable to use this source term");
    }
    auto const& blank =
        this->m_sim.repo().get_field("eb_blank")(lev).const_arrays();

    const auto& geom = m_mesh.Geom(lev);
    const auto& problo = geom.ProbLoArray();
    const auto& dx = geom.CellSizeArray();
    const amrex::Real drag_coefficient = m_drag_coefficient;
    const amrex::Real Cd = drag_coefficient / dx[2];
    const amrex::Real cd_max = 10.0_rt;
    const amrex::Real T0_const = m_soil_temperature;

    // Use IDW model if enabled
    if (m_use_idw_model) {
        const amrex::Real* idw_x_ptr = m_idw_x.data();
        const amrex::Real* idw_y_ptr = m_idw_y.data();
        const amrex::Real* idw_z_ptr = m_idw_z.data();
        const amrex::Real* idw_temp_ptr = m_idw_temp.data();
        const int num_points = m_idw_num_points;

        amrex::ParallelFor(
            src_term, amrex::IntVect(0), AMREX_SPACEDIM,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int /*n*/) {
                const amrex::Real ux1 = vel[nbx](i, j, k, 0);
                const amrex::Real uy1 = vel[nbx](i, j, k, 1);
                const amrex::Real uz1 = vel[nbx](i, j, k, 2);
                const amrex::Real theta = temperature[nbx](i, j, k, 0);
                const amrex::Real m =
                    std::sqrt((ux1 * ux1) + (uy1 * uy1) + (uz1 * uz1));
                const amrex::Real CdM = std::min(
                    Cd / (m + kynema_sgf::constants::EPS), cd_max / dx[2]);

                // Compute cell center coordinates
                const amrex::Real x = problo[0] + (i + 0.5_rt) * dx[0];
                const amrex::Real y = problo[1] + (j + 0.5_rt) * dx[1];
                const amrex::Real z = problo[2] + (k + 0.5_rt) * dx[2];

                // Compute T0 using inverse distance weighting
                amrex::Real sum_weights = 0.0_rt;
                amrex::Real weighted_temp = 0.0_rt;

                for (int ip = 0; ip < num_points; ++ip) {
                    const amrex::Real dist_x = x - idw_x_ptr[ip];
                    const amrex::Real dist_y = y - idw_y_ptr[ip];
                    const amrex::Real dist_z = z - idw_z_ptr[ip];
                    const amrex::Real dist = std::sqrt(
                        dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);

                    // Avoid division by zero for points very close to IDW
                    // points
                    const amrex::Real weight =
                        1.0_rt / (dist + kynema_sgf::constants::EPS);
                    sum_weights += weight;
                    weighted_temp += weight * idw_temp_ptr[ip];
                }

                const amrex::Real T0 = weighted_temp / sum_weights;

                src_arrs[nbx](i, j, k, 0) -=
                    (CdM * (theta - T0) * blank[nbx](i, j, k, 0));
            });
    } else if (m_use_heatflux_model) {
        // Use heat flux based soil temperature model
        // T0 is computed from: Q = k * (T_surface - T0) / h
        // Therefore: T0 = T_surface - (Q * h / k)

        const amrex::Real heatflux = m_surface_heatflux;
        const amrex::Real conductivity = m_thermal_conductivity;
        const amrex::Real thickness = m_boundary_layer_height;

        amrex::ParallelFor(
            src_term, amrex::IntVect(0), AMREX_SPACEDIM,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int /*n*/) {
                const amrex::Real ux1 = vel[nbx](i, j, k, 0);
                const amrex::Real uy1 = vel[nbx](i, j, k, 1);
                const amrex::Real uz1 = vel[nbx](i, j, k, 2);
                const amrex::Real theta = temperature[nbx](i, j, k, 0);
                const amrex::Real m =
                    std::sqrt((ux1 * ux1) + (uy1 * uy1) + (uz1 * uz1));
                const amrex::Real CdM = std::min(
                    Cd / (m + kynema_sgf::constants::EPS), cd_max / dx[2]);

                // Compute T0 from heat flux
                // Assumes theta is the surface temperature at this location
                const amrex::Real T0 =
                    theta - (heatflux * thickness / conductivity);

                src_arrs[nbx](i, j, k, 0) -=
                    (CdM * (theta - T0) * blank[nbx](i, j, k, 0));
            });
    } else {
        // Use constant soil temperature model
        amrex::ParallelFor(
            src_term, amrex::IntVect(0), AMREX_SPACEDIM,
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int /*n*/) {
                const amrex::Real ux1 = vel[nbx](i, j, k, 0);
                const amrex::Real uy1 = vel[nbx](i, j, k, 1);
                const amrex::Real uz1 = vel[nbx](i, j, k, 2);
                const amrex::Real theta = temperature[nbx](i, j, k, 0);
                const amrex::Real m =
                    std::sqrt((ux1 * ux1) + (uy1 * uy1) + (uz1 * uz1));
                const amrex::Real CdM = std::min(
                    Cd / (m + kynema_sgf::constants::EPS), cd_max / dx[2]);
                src_arrs[nbx](i, j, k, 0) -=
                    (CdM * (theta - T0_const) * blank[nbx](i, j, k, 0));
            });
    }
}

} // namespace kynema_sgf::pde::temperature
