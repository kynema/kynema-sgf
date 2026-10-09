#include "src/transport_models/TransportModel.H"
#include "src/wind_energy/ABLWallFunction.H"
#include "src/wind_energy/ABL.H"
#include "src/utilities/tensor_ops.H"
#include "src/utilities/trig_ops.H"
#include "src/diffusion/diffusion.H"
#include "src/wind_energy/ShearStress.H"
#include "src/wind_energy/MOData.H"
#include "src/utilities/linear_interpolation.H"

#include <cmath>
#include <type_traits>

#include "AMReX_ParmParse.H"
#include "AMReX_Print.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf {

ABLWallFunction::ABLWallFunction(const CFDSim& sim)
    : m_sim(sim), m_mesh(sim.mesh())
{
    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }

    amrex::ParmParse pp("ABL");
    pp.query("kappa", m_mo.kappa);
    pp.query("mo_gamma_m", m_mo.gamma_m);
    pp.query("mo_gamma_h", m_mo.gamma_h);
    pp.query("mo_beta_m", m_mo.beta_m);
    pp.query("mo_beta_h", m_mo.beta_h);
    const char* z0_same = "surface_roughness_z0";
    const char* z0_aero = "aerodynamic_roughness_length";
    const char* z0_therm = "thermal_roughness_length";
    pp.query(z0_same, m_mo.z0);
    if (pp.contains(z0_aero) && pp.contains(z0_therm)) {
        pp.query(z0_aero, m_mo.z0);
        pp.query(z0_therm, m_mo.z0t);
        if (pp.contains(z0_same)) {
            amrex::Print()
                << "WARNING: ABLWallFunction: Both " << z0_same
                << " and separate roughness lengths (aerodynamic "
                   "and thermal) are specified. Using separate roughness "
                   "lengths and ignoring combined parameter "
                << z0_same << ".\n";
        }
    } else if (
        pp.contains(z0_same) &&
        (pp.contains(z0_aero) || pp.contains(z0_therm))) {
        amrex::Abort(
            "ABLWallFunction parameter conflict: Roughness lengths must be "
            "specified as the same (" +
            std::string(z0_same) + ") or as different (" +
            std::string(z0_aero) + " and " + std::string(z0_therm) + ").");
    } else {
        m_mo.z0t = m_mo.z0;
    }
    pp.query("normal_direction", m_direction);
    AMREX_ASSERT((0 <= m_direction) && (m_direction < AMREX_SPACEDIM));

    if (pp.contains("log_law_height")) {
        m_use_fch = false;
        pp.get("log_law_height", m_mo.zref);
    } else {
        m_use_fch = true;
        amrex::Print()
            << "ABLWallFunction: log_law_height not specified for ABL physics. "
               "Assuming log_law_height = first cell height"
            << '\n';
    }
    m_wall_pos = m_sim.mesh().Geom(0).ProbLo(m_direction);
    pp.query("wall_position", m_wall_pos);
    pp.query("level_means_method", m_level_means_method);
    m_level_means_method = amrex::toLower(m_level_means_method);
    if ((m_level_means_method != "mo_profile") &&
        (m_level_means_method != "none")) {
        amrex::Abort(
            "ABL.level_means_method must be mo_profile or none, got " +
            m_level_means_method);
    }

    if (pp.contains("surface_temp_flux")) {
        pp.query("surface_temp_flux", m_mo.surf_temp_flux);
        amrex::Print()
            << "ABLWallFunction: Surface temperature flux mode is selected."
            << '\n';
    } else if (pp.contains("surface_temp_timetable")) {
        pp.query("surface_temp_timetable", m_surf_temp_timetable);
        m_tempflux = false;
        m_temp_table = true;
        amrex::Print() << "ABLWallFunction: Surface temperature time table "
                          "mode is selected."
                       << '\n';
        if (!m_surf_temp_timetable.empty()) {
            std::ifstream ifh(m_surf_temp_timetable, std::ios::in);
            if (!ifh.good()) {
                amrex::Abort(
                    "Cannot find surface_temp_timetable file: " +
                    m_surf_temp_timetable);
            }
            amrex::Real data_time;
            amrex::Real data_value;
            ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            while (ifh >> data_time) {
                ifh >> data_value;
                m_surf_temp_time.push_back(data_time);
                m_surf_temp_value.push_back(data_value);
            }
        }
    } else if (pp.contains("surface_temp_rate")) {
        m_tempflux = false;
        pp.get("surface_temp_rate", m_surf_temp_rate);
        amrex::Print()
            << "ABLWallFunction: Surface temperature rate mode is selected."
            << '\n';
        if (pp.contains("surface_temp_init")) {
            pp.get("surface_temp_init", m_surf_temp_init);
        } else {
            amrex::Print()
                << "ABLWallFunction: Initial surface temperature not found for "
                   "ABL. Assuming to be equal to the reference temperature"
                << '\n';
            m_surf_temp_init = sim.transport_model().reference_temperature();
        }
        if (pp.contains("surface_temp_rate_tstart")) {
            pp.get("surface_temp_rate_tstart", m_surf_temp_rate_tstart);
        } else {
            amrex::Print()
                << "ABLWallFunction: Surface temperature heating/cooling start "
                   "time (surface_temp_rate_tstart) not found for ABL. "
                   "Assuming zero."
                << m_surf_temp_rate_tstart << '\n';
        }

    } else {
        amrex::Print() << "ABLWallFunction: Neither surface_temp_flux nor "
                          "surface_temp_rate specified for ABL physics. "
                          "Assuming Neutral Stratification"
                       << '\n';
    }

    if (pp.contains("inflow_outflow_mode")) {
        pp.query("inflow_outflow_mode", m_inflow_outflow);
        if (m_inflow_outflow) {
            pp.query("wf_velocity", m_wf_vel);
            pp.query("wf_vmag", m_wf_vmag);
            pp.query("wf_theta", m_wf_theta);
            amrex::Print() << "ABLWallFunction: Inflow/Outflow mode is turned "
                              "on. Please make sure wall shear stress type is "
                              "set to local."
                           << '\n';
        }
    }

    m_mo.alg_type = m_tempflux ? MOData::ThetaCalcType::HEAT_FLUX
                               : MOData::ThetaCalcType::SURFACE_TEMPERATURE;
    m_mo.gravity = utils::vec_mag(m_gravity.data());
}

void ABLWallFunction::init_log_law_height(int max_level)
{
    if (m_use_fch) {
        const auto& geom = m_mesh.Geom(max_level);
        m_mo.zref = 0.5_rt * geom.CellSize(m_direction);
    }
}

void ABLWallFunction::update_umean(
    const VelPlaneAveraging& vpa, const FieldPlaneAveraging& tpa)
{
    const auto& time = m_sim.time();

    if (!m_tempflux) {
        if (!m_temp_table) {
            m_mo.surf_temp =
                m_surf_temp_init +
                (m_surf_temp_rate *
                 amrex::max<amrex::Real>(
                     time.current_time() - m_surf_temp_rate_tstart, 0.0_rt) /
                 3600.0_rt);
        } else {
            m_mo.surf_temp = kynema_sgf::interp::linear(
                m_surf_temp_time, m_surf_temp_value, time.current_time());
        }
        amrex::Print() << "Current surface temperature: " << m_mo.surf_temp
                       << '\n';
    }

    if (m_inflow_outflow) {
        m_mo.vel_mean[0] = m_wf_vel[0];
        m_mo.vel_mean[1] = m_wf_vel[1];
        m_mo.vmag_mean = m_wf_vmag;
        m_mo.Su_mean = 0.0_rt; // TODO: need to fill this correctly
        m_mo.Sv_mean = 0.0_rt; // TODO: need to fill this correctly
        m_mo.theta_mean = m_wf_theta;
    } else {
        m_mo.vel_mean[0] =
            vpa.line_average_interpolated(m_wall_pos + m_mo.zref, 0);
        m_mo.vel_mean[1] =
            vpa.line_average_interpolated(m_wall_pos + m_mo.zref, 1);
        m_mo.vmag_mean =
            vpa.line_hvelmag_average_interpolated(m_wall_pos + m_mo.zref);
        m_mo.Su_mean = vpa.line_su_average_interpolated(m_wall_pos + m_mo.zref);
        m_mo.Sv_mean = vpa.line_sv_average_interpolated(m_wall_pos + m_mo.zref);
        m_mo.theta_mean =
            tpa.line_average_interpolated(m_wall_pos + m_mo.zref, 0);
    }

    m_mo.update_fluxes();

    // Mean quantities of every mesh level that touches the wall
    if (m_inflow_outflow) {
        m_mo_lev.clear();
    } else {
        update_level_means();
    }
}

void ABLWallFunction::update_level_means()
{
    m_mo_lev.clear();
    if (m_level_means_method == "none") {
        return;
    }

    const auto& repo = m_sim.repo();
    const int nlevels = repo.num_active_levels();
    // A single-level mesh keeps the plane averages at the reference height
    if (nlevels < 2) {
        return;
    }

    // The means only serve the wall models at zlo
    const amrex::Orientation zlo(amrex::Direction::z, amrex::Orientation::low);
    const bool vel_wall =
        repo.get_field("velocity").bc_type()[zlo] == BC::wall_model;
    const bool temp_wall =
        repo.field_exists("temperature") &&
        (repo.get_field("temperature").bc_type()[zlo] == BC::wall_model);
    if (!vel_wall && !temp_wall) {
        return;
    }

    // Same wall as ABLVelWallFunc and ABLTempWallFunc, which apply the wall
    // models at zlo. Layer of the wall-adjacent cells of a level:
    constexpr int idim = 2;
    const auto wall_layer = [this](const int lev) {
        amrex::Box bx = m_mesh.Geom(lev).Domain();
        bx.setBig(idim, bx.smallEnd(idim));
        return bx;
    };

    // Refinement that does not reach the wall leaves level 0 alone on the
    // wall, which then keeps the plane averages at the reference height as
    // on a single-level mesh. The box arrays tell without touching any data
    bool fine_level_on_wall = false;
    for (int lev = 1; lev < nlevels; ++lev) {
        fine_level_on_wall = fine_level_on_wall ||
                             m_mesh.boxArray(lev).intersects(wall_layer(lev));
    }
    if (!fine_level_on_wall) {
        return;
    }

    // The wall models combine the local values of the first cell of each
    // level with the plane averages entering the Monin-Obukhov data. On a
    // mesh where a finer level touches the wall, the first cells of the
    // levels sit at different heights, so each level gets the plane averages
    // of the reference height carried to its own first-cell height z_l
    // along the Monin-Obukhov profile of the surface layer:
    //
    //   <u>_l     = <u>     phi_m(z_l) / phi_m(z_ref)      (also v, |u_h|)
    //   <theta>_l = theta_s + (<theta> - theta_s) phi_h(z_l) / phi_h(z_ref)
    //
    // with phi_m(z) = ln(z/z0) - psi_m(z/L) and phi_h(z) = ln(z/z0t) -
    // psi_h(z/L). When the plane-mean profiles follow the Monin-Obukhov
    // profile, the mean stress and heat flux of every level then equal those
    // of the reference height. The friction velocity, the Obukhov length and
    // the surface heat flux stay those of the reference height.
    const amrex::Real phi_m_ref = m_mo.phi_m();
    const amrex::Real phi_h_ref = m_mo.phi_h();
    m_mo_lev.resize(nlevels, m_mo);
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = m_mesh.Geom(lev);
        if (!m_mesh.boxArray(lev).intersects(wall_layer(lev))) {
            // No wall cells on this level: keeps the copy of m_mo
            continue;
        }

        // Height of the first cell center above the wall, measured from the
        // wall position like the reference height
        const amrex::Real zref_lev =
            geom.ProbLo(idim) + (0.5_rt * geom.CellSize(idim)) - m_wall_pos;
        // The log law does not hold at or below the roughness height, and
        // in strongly unstable conditions phi_m or phi_h can vanish or change
        // sign just above it: such a level keeps the means at the reference
        // height, as without refinement
        const bool above_roughness = zref_lev > amrex::max(m_mo.z0, m_mo.z0t);
        const amrex::Real phi_m_lev =
            above_roughness ? m_mo.phi_m(zref_lev) : 0.0_rt;
        const amrex::Real phi_h_lev =
            above_roughness ? m_mo.phi_h(zref_lev) : 0.0_rt;
        if (!(phi_m_lev > 0.0_rt) || !(phi_h_lev > 0.0_rt)) {
            if (!m_warned_low_first_cell) {
                amrex::Print()
                    << "WARNING: ABLWallFunction: the first cell of level "
                    << lev << " sits at " << zref_lev
                    << " m, where the log-law profile functions are not "
                       "positive; this level uses the plane averages at the "
                       "reference height\n";
                m_warned_low_first_cell = true;
            }
            continue;
        }

        // The surface temperature, friction velocity, Obukhov length and
        // surface heat flux are those of the reference height. With the
        // scaled mean temperature, kappa u* (theta_s - <theta>_l) / phi_h(z_l)
        // is the reference-height heat flux on every level
        auto& mo_lev = m_mo_lev[lev];
        mo_lev.zref = zref_lev;
        const amrex::Real rm = phi_m_lev / phi_m_ref;
        const amrex::Real rh = phi_h_lev / phi_h_ref;
        mo_lev.vel_mean[0] = m_mo.vel_mean[0] * rm;
        mo_lev.vel_mean[1] = m_mo.vel_mean[1] * rm;
        mo_lev.vmag_mean = m_mo.vmag_mean * rm;
        mo_lev.Su_mean = m_mo.Su_mean * rm * rm;
        mo_lev.Sv_mean = m_mo.Sv_mean * rm * rm;
        mo_lev.theta_mean =
            m_mo.surf_temp + ((m_mo.theta_mean - m_mo.surf_temp) * rh);
    }
}

void ABLWallFunction::update_tflux(const amrex::Real tflux)
{
    m_mo.surf_temp_flux = tflux;
    // The per-level data carry the same surface heat flux. Their surface
    // temperatures, like that of m_mo, stay those of the last update_umean
    for (auto& mo_lev : m_mo_lev) {
        mo_lev.surf_temp_flux = tflux;
    }
}

ABLVelWallFunc::ABLVelWallFunc(
    Field& /*unused*/, const ABLWallFunction& wall_func)
    : m_wall_func(wall_func)
{
    amrex::ParmParse pp("ABL");
    pp.query("wall_shear_stress_type", m_wall_shear_stress_type);
    pp.query("wall_het_model", m_wall_het_model);
    pp.query("monin_obukhov_length", m_monin_obukhov_length);
    m_wall_shear_stress_type = amrex::toLower(m_wall_shear_stress_type);

    if (m_wall_shear_stress_type == "constant" ||
        m_wall_shear_stress_type == "local" ||
        m_wall_shear_stress_type == "schumann" ||
        m_wall_shear_stress_type == "donelan" ||
        m_wall_shear_stress_type == "moeng") {
        amrex::Print() << "Shear Stress model: " << m_wall_shear_stress_type
                       << '\n';
    } else {
        amrex::Abort("Shear Stress wall model input mistake");
    }
}

namespace {
//! Whether a wall model takes the mean quantities of each mesh level. The
//! Donelan model selects its drag coefficient (and, with a specified
//! surface temperature, its heat flux) from the mean wind at the reference
//! height, ABL.log_law_height, so it keeps the reference-height data on
//! every level
template <typename WallModel>
constexpr bool uses_level_means()
{
    return !std::is_same_v<WallModel, ShearStressDonelan>;
}
} // namespace

template <typename ShearStress>
void ABLVelWallFunc::wall_model(Field& velocity, const FieldState rho_state)
{
    BL_PROFILE("kynema-sgf::ABLVelWallFunc");

    constexpr int idim = 2;
    const auto& repo = velocity.repo();
    const auto& density = repo.get_field("density", rho_state);
    const auto& viscosity = repo.get_field("velocity_mueff");
    const int nlevels = repo.num_active_levels();
    amrex::Orientation zlo(amrex::Direction::z, amrex::Orientation::low);
    amrex::Orientation zhi(amrex::Direction::z, amrex::Orientation::high);
    if (velocity.bc_type()[zhi] == BC::wall_model) {
        amrex::Abort("ABL wall models are not applicable to a zhi BC");
    }
    if (velocity.bc_type()[zlo] != BC::wall_model) {
        return;
    }
    const auto& mo = m_wall_func.mo();
    const bool has_terrain = repo.int_field_exists("terrain_blank");
    const auto* m_terrain_blank =
        has_terrain ? &repo.get_int_field("terrain_blank") : nullptr;
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = repo.mesh().Geom(lev);
        const auto& domain = geom.Domain();
        amrex::MFItInfo mfi_info{};
        const auto& dx = geom.CellSizeArray();
        const auto& rho_lev = density(lev);
        const auto& vold_lev = velocity.state(FieldState::Old)(lev);
        auto& vel_lev = velocity(lev);
        const auto& eta_lev = viscosity(lev);
        // Shear-stress model with the mean quantities of this level, or of
        // the reference height for the Donelan model
        const ShearStress tau(
            uses_level_means<ShearStress>() ? m_wall_func.mo(lev)
                                            : m_wall_func.mo());

        if (amrex::Gpu::notInLaunchRegion()) {
            mfi_info.SetDynamic(true);
        }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(vel_lev, mfi_info); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.validbox();
            const auto& varr = vel_lev.array(mfi);
            const auto& vold_arr = vold_lev.const_array(mfi);
            const auto& den = rho_lev.const_array(mfi);
            const auto& eta = eta_lev.const_array(mfi);
            const auto& blank_arr =
                has_terrain ? (*m_terrain_blank)(lev).const_array(mfi)
                            : amrex::Array4<int>();
            if (bx.smallEnd(idim) == domain.smallEnd(idim) &&
                velocity.bc_type()[zlo] == BC::wall_model) {
                if (m_wall_het_model == "mol") {
                    const amrex::Real z = 0.5_rt * dx[2];
                    const amrex::Real zeta = z / m_monin_obukhov_length;
                    const amrex::Real psi_m = mo.calc_psi_m(zeta);
                    const amrex::Real z0 = mo.z0;
                    const amrex::Real kappa = mo.kappa;
                    amrex::ParallelFor(
                        amrex::bdryLo(bx, idim),
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                            const amrex::Real mu = eta(i, j, k);
                            const amrex::Real uu = vold_arr(i, j, k, 0);
                            const amrex::Real vv = vold_arr(i, j, k, 1);
                            const amrex::Real wspd =
                                std::sqrt((uu * uu) + (vv * vv));
                            const amrex::Real drag = std::log(z / z0) - psi_m;
                            const amrex::Real ustar = wspd * kappa / drag;
                            // Dirichlet BC
                            varr(i, j, k - 1, 2) = 0.0_rt;
                            const amrex::Real blankTerrain =
                                (has_terrain) ? 1 - blank_arr(i, j, k, 0)
                                              : 1.0_rt;
                            // Shear stress BC
                            // Blank Terrain added to keep the boundary
                            // condition backward compatible while adding
                            // terrain sensitive BC
                            varr(i, j, k - 1, 0) = blankTerrain * ustar *
                                                   ustar * uu / wspd *
                                                   den(i, j, k) / mu;
                            varr(i, j, k - 1, 1) = blankTerrain * ustar *
                                                   ustar * vv / wspd *
                                                   den(i, j, k) / mu;
                        });
                } else {
                    amrex::ParallelFor(
                        amrex::bdryLo(bx, idim),
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                            const amrex::Real mu = eta(i, j, k);
                            const amrex::Real uu = vold_arr(i, j, k, 0);
                            const amrex::Real vv = vold_arr(i, j, k, 1);
                            const amrex::Real wspd =
                                std::sqrt((uu * uu) + (vv * vv));

                            // Dirichlet BC
                            varr(i, j, k - 1, 2) = 0.0_rt;
                            const amrex::Real blankTerrain =
                                (has_terrain) ? 1 - blank_arr(i, j, k, 0)
                                              : 1.0_rt;
                            // Shear stress BC
                            // Blank Terrain added to keep the boundary
                            // condition backward compatible while adding
                            // terrain sensitive BC
                            varr(i, j, k - 1, 0) = blankTerrain *
                                                   tau.calc_vel_x(uu, wspd) *
                                                   den(i, j, k) / mu;
                            varr(i, j, k - 1, 1) = blankTerrain *
                                                   tau.calc_vel_y(vv, wspd) *
                                                   den(i, j, k) / mu;
                        });
                }
            }
        }
    }
}

void ABLVelWallFunc::operator()(Field& velocity, const FieldState rho_state)
{
    if (m_wall_shear_stress_type == "moeng") {
        wall_model<ShearStressMoeng>(velocity, rho_state);
    } else if (m_wall_shear_stress_type == "constant") {
        wall_model<ShearStressConstant>(velocity, rho_state);
    } else if (m_wall_shear_stress_type == "local") {
        wall_model<ShearStressLocal>(velocity, rho_state);
    } else if (m_wall_shear_stress_type == "schumann") {
        wall_model<ShearStressSchumann>(velocity, rho_state);
    } else if (m_wall_shear_stress_type == "donelan") {
        wall_model<ShearStressDonelan>(velocity, rho_state);
    }
}

ABLTempWallFunc::ABLTempWallFunc(
    Field& /*unused*/, const ABLWallFunction& wall_fuc)
    : m_wall_func(wall_fuc)
{
    amrex::ParmParse pp("ABL");
    pp.query("wall_shear_stress_type", m_wall_shear_stress_type);
    pp.query("wall_het_model", m_wall_het_model);
    pp.query("monin_obukhov_length", m_monin_obukhov_length);
    m_wall_shear_stress_type = amrex::toLower(m_wall_shear_stress_type);
    amrex::Print() << "Heat Flux model: " << m_wall_shear_stress_type << '\n';
}

template <typename HeatFlux>
void ABLTempWallFunc::wall_model(Field& temperature, const FieldState rho_state)
{
    constexpr int idim = 2;
    auto& repo = temperature.repo();
    const auto& mo = m_wall_func.mo();

    // Return early if the user hasn't requested a wall model BC for temperature
    amrex::Orientation zlo(amrex::Direction::z, amrex::Orientation::low);
    amrex::Orientation zhi(amrex::Direction::z, amrex::Orientation::high);
    if (temperature.bc_type()[zhi] == BC::wall_model) {
        amrex::Abort("ABL wall models are not applicable to a zhi BC");
    }
    if (temperature.bc_type()[zlo] != BC::wall_model) {
        return;
    }

    BL_PROFILE("kynema-sgf::ABLTempWallFunc");
    auto& velocity = repo.get_field("velocity");
    const auto& density = repo.get_field("density", rho_state);
    const auto& alpha = repo.get_field("temperature_mueff");
    const int nlevels = repo.num_active_levels();
    const bool has_terrain = repo.int_field_exists("terrain_blank");
    const auto* m_terrain_blank =
        has_terrain ? &repo.get_int_field("terrain_blank") : nullptr;
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = repo.mesh().Geom(lev);
        const auto& domain = geom.Domain();
        amrex::MFItInfo mfi_info{};
        const auto& dx = geom.CellSizeArray();
        const auto& rho_lev = density(lev);
        const auto& vold_lev = velocity.state(FieldState::Old)(lev);
        const auto& told_lev = temperature.state(FieldState::Old)(lev);
        auto& theta = temperature(lev);
        const auto& eta_lev = alpha(lev);
        // Heat-flux model with the mean quantities of this level, or of the
        // reference height for the Donelan model
        const HeatFlux tau(
            uses_level_means<HeatFlux>() ? m_wall_func.mo(lev)
                                         : m_wall_func.mo());

        if (amrex::Gpu::notInLaunchRegion()) {
            mfi_info.SetDynamic(true);
        }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(theta, mfi_info); mfi.isValid(); ++mfi) {
            const auto& bx = mfi.validbox();
            const auto& vold_arr = vold_lev.const_array(mfi);
            const auto& told_arr = told_lev.const_array(mfi);
            const auto& tarr = theta.array(mfi);
            const auto& den = rho_lev.const_array(mfi);
            const auto& eta = eta_lev.const_array(mfi);
            const auto& blank_arr =
                has_terrain ? (*m_terrain_blank)(lev).const_array(mfi)
                            : amrex::Array4<int>();

            if (bx.smallEnd(idim) == domain.smallEnd(idim) &&
                temperature.bc_type()[zlo] == BC::wall_model) {
                if (m_wall_het_model == "mol") {
                    const amrex::Real z = 0.5_rt * dx[2];
                    const amrex::Real zeta = z / m_monin_obukhov_length;
                    const amrex::Real psi_m = mo.calc_psi_m(zeta);
                    const amrex::Real z0 = mo.z0;
                    const amrex::Real kappa = mo.kappa;
                    const amrex::Real gravity_mod = 9.81_rt;
                    const amrex::Real monin_obukhov_length =
                        m_monin_obukhov_length;
                    amrex::ParallelFor(
                        amrex::bdryLo(bx, idim),
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                            const amrex::Real alphaT = eta(i, j, k);
                            const amrex::Real uu = vold_arr(i, j, k, 0);
                            const amrex::Real vv = vold_arr(i, j, k, 1);
                            const amrex::Real wspd =
                                std::sqrt((uu * uu) + (vv * vv));
                            const amrex::Real theta2 = told_arr(i, j, k);
                            const amrex::Real drag = std::log(z / z0) - psi_m;
                            const amrex::Real ustar = wspd * kappa / drag;
                            const amrex::Real thetastar =
                                theta2 * ustar * ustar /
                                (kappa * gravity_mod * monin_obukhov_length);
                            const amrex::Real surf_temp_flux =
                                ustar * thetastar;
                            const amrex::Real blankTerrain =
                                (has_terrain) ? 1 - blank_arr(i, j, k, 0)
                                              : 1.0_rt;
                            tarr(i, j, k - 1) = blankTerrain * den(i, j, k) *
                                                surf_temp_flux / alphaT;
                        });
                } else {
                    amrex::ParallelFor(
                        amrex::bdryLo(bx, idim),
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                            const amrex::Real alphaT = eta(i, j, k);
                            const amrex::Real uu = vold_arr(i, j, k, 0);
                            const amrex::Real vv = vold_arr(i, j, k, 1);
                            const amrex::Real wspd =
                                std::sqrt((uu * uu) + (vv * vv));
                            const amrex::Real theta2 = told_arr(i, j, k);
                            const amrex::Real blankTerrain =
                                (has_terrain) ? 1 - blank_arr(i, j, k, 0)
                                              : 1.0_rt;
                            tarr(i, j, k - 1) = blankTerrain * den(i, j, k) *
                                                tau.calc_theta(wspd, theta2) /
                                                alphaT;
                        });
                }
            }
        }
    }
}

void ABLTempWallFunc::operator()(Field& temperature, const FieldState rho_state)
{
    if (m_wall_shear_stress_type == "moeng") {
        wall_model<ShearStressMoeng>(temperature, rho_state);
    } else if (m_wall_shear_stress_type == "constant") {
        wall_model<ShearStressConstant>(temperature, rho_state);
    } else if (m_wall_shear_stress_type == "local") {
        wall_model<ShearStressLocal>(temperature, rho_state);
    } else if (m_wall_shear_stress_type == "schumann") {
        wall_model<ShearStressSchumann>(temperature, rho_state);
    } else if (m_wall_shear_stress_type == "donelan") {
        wall_model<ShearStressDonelan>(temperature, rho_state);
    }
}

} // namespace kynema_sgf
