#include <cmath>
#include <limits>

#include "src/turbulence/RANS/KLAxellSeparation.H"
#include "src/equation_systems/PDEBase.H"
#include "src/turbulence/TurbModelDefs.H"
#include "src/fvm/gradient.H"
#include "src/fvm/strainrate.H"
#include "src/utilities/math_ops.H"
#include "AMReX_ParmParse.H"

using namespace amrex::literals;

namespace kynema_sgf {
namespace turbulence {

template <typename Transport>
KLAxellSeparation<Transport>::KLAxellSeparation(CFDSim& sim)
    : KLAxell<Transport>(sim)
{
    amrex::ParmParse pp("KLAxellSeparation");
    pp.query("pressure_gradient_sensor", m_use_pressure_gradient_sensor);
    pp.query("realizable_cmu", m_use_realizable_cmu);
    pp.query("sensor_source", m_sensor_source);
    if (m_sensor_source != "pressure" && m_sensor_source != "velocity") {
        amrex::Abort(
            "KLAxellSeparation.sensor_source must be pressure or velocity");
    }
    if (m_use_realizable_cmu && !m_use_pressure_gradient_sensor) {
        amrex::Abort(
            "KLAxellSeparation.realizable_cmu requires "
            "KLAxellSeparation.pressure_gradient_sensor = true");
    }
    // The production cap is applied by the TKE source (KransAxell)
    bool use_production_cap = false;
    pp.query("production_cap", use_production_cap);
    if (use_production_cap && !m_use_pressure_gradient_sensor) {
        amrex::Abort(
            "KLAxellSeparation.production_cap requires "
            "KLAxellSeparation.pressure_gradient_sensor = true");
    }
    // The destruction boost is applied by the TKE source (KransAxell)
    bool use_destruction_boost = false;
    pp.query("destruction_boost", use_destruction_boost);
    if (use_destruction_boost && !m_use_pressure_gradient_sensor) {
        amrex::Abort(
            "KLAxellSeparation.destruction_boost requires "
            "KLAxellSeparation.pressure_gradient_sensor = true");
    }
    pp.query("curvature_correction", m_use_curvature_correction);
    pp.query("curvature_model", m_curvature_model);
    if (m_curvature_model != "rotation_function" &&
        m_curvature_model != "richardson") {
        amrex::Abort(
            "KLAxellSeparation.curvature_model must be rotation_function or "
            "richardson");
    }
    if (m_use_pressure_gradient_sensor) {
        m_pressure_gradient_sensor =
            &sim.repo().declare_field("pressure_gradient_sensor", 1);
    }

    // The relaxation time decides whether the gate field exists, so it is
    // read here rather than in parse_model_coeffs. The gate needs the sensor:
    // without it the default relaxation time has no effect, and a positive
    // value set in the input is an error.
    amrex::ParmParse pp_coeffs("KLAxellSeparation_coeffs");
    const bool relaxation_time_set =
        pp_coeffs.query("gate_relaxation_time", m_gate_relaxation_time) != 0;
    if (m_gate_relaxation_time < 0.0_rt) {
        amrex::Abort(
            "KLAxellSeparation_coeffs.gate_relaxation_time must not be "
            "negative");
    }
    if (relaxation_time_set && (m_gate_relaxation_time > 0.0_rt) &&
        !m_use_pressure_gradient_sensor) {
        amrex::Abort(
            "KLAxellSeparation_coeffs.gate_relaxation_time requires "
            "KLAxellSeparation.pressure_gradient_sensor = true");
    }
    if (m_use_pressure_gradient_sensor && (m_gate_relaxation_time > 0.0_rt)) {
        m_separation_gate = &sim.repo().declare_field("separation_gate", 1, 1);
        m_separation_gate->set_default_fillpatch_bc(sim.time());
        m_separation_gate->fillpatch_on_regrid() = true;
        sim.io_manager().register_restart_var("separation_gate");
    }
}

template <typename Transport>
KLAxellSeparation<Transport>::~KLAxellSeparation() = default;

template <typename Transport>
void KLAxellSeparation<Transport>::parse_model_coeffs()
{
    KLAxell<Transport>::parse_model_coeffs();
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("sensor_velocity_weight", m_sensor_velocity_weight);
    pp.query("sensor_threshold", m_sensor_threshold);
    pp.query("realizable_cmu_strength", m_realizable_cmu_strength);
    pp.query("production_cap_ratio", m_production_cap_ratio);
    pp.query("destruction_boost_factor", m_destruction_boost_factor);
    pp.query("curvature_coefficient", m_curvature_coefficient);
    pp.query("curvature_cr1", m_curvature_cr1);
    pp.query("curvature_cr2", m_curvature_cr2);
    pp.query("curvature_cr3", m_curvature_cr3);
    if (m_curvature_coefficient < 0.0_rt) {
        amrex::Abort(
            "KLAxellSeparation_coeffs.curvature_coefficient must not be "
            "negative");
    }
    // A negative weight or strength could make the sensor scale or the
    // limiter denominator vanish or change sign
    if (m_sensor_velocity_weight < 0.0_rt) {
        amrex::Abort(
            "KLAxellSeparation_coeffs.sensor_velocity_weight must not be "
            "negative");
    }
    if (m_realizable_cmu_strength < 0.0_rt) {
        amrex::Abort(
            "KLAxellSeparation_coeffs.realizable_cmu_strength must not be "
            "negative");
    }
    if (m_sensor_threshold <= 0.0_rt) {
        amrex::Abort(
            "KLAxellSeparation_coeffs.sensor_threshold must be positive");
    }
    if (m_production_cap_ratio <= 0.0_rt) {
        amrex::Abort(
            "KLAxellSeparation_coeffs.production_cap_ratio must be positive");
    }
    if (m_destruction_boost_factor < 1.0_rt) {
        amrex::Abort(
            "KLAxellSeparation_coeffs.destruction_boost_factor must not be "
            "below 1");
    }
}

template <typename Transport>
TurbulenceModel::CoeffsDictType
KLAxellSeparation<Transport>::model_coeffs() const
{
    auto coeffs = KLAxell<Transport>::model_coeffs();
    coeffs["sensor_velocity_weight"] = m_sensor_velocity_weight;
    coeffs["sensor_threshold"] = m_sensor_threshold;
    coeffs["realizable_cmu_strength"] = m_realizable_cmu_strength;
    coeffs["gate_relaxation_time"] = m_gate_relaxation_time;
    coeffs["production_cap_ratio"] = m_production_cap_ratio;
    coeffs["destruction_boost_factor"] = m_destruction_boost_factor;
    coeffs["curvature_coefficient"] = m_curvature_coefficient;
    coeffs["curvature_cr1"] = m_curvature_cr1;
    coeffs["curvature_cr2"] = m_curvature_cr2;
    coeffs["curvature_cr3"] = m_curvature_cr3;
    return coeffs;
}

template <typename Transport>
void KLAxellSeparation<Transport>::post_init_actions()
{
    KLAxell<Transport>::post_init_actions();
    m_geometry = this->m_sim.repo().create_scratch_field(
        klaxell_separation::geom_ncomp, 1);
    // The curvature Rt is read by update_alphaeff before the first viscosity
    // update
    for (int lev = 0; lev < this->m_sim.repo().num_active_levels(); ++lev) {
        (*m_geometry)(lev).setVal(0.0_rt);
    }
    if (m_use_pressure_gradient_sensor) {
        m_pressure_gradient_sensor->setVal(0.0_rt);
    }
    if (m_use_realizable_cmu) {
        m_strain = this->m_sim.repo().create_scratch_field(1, 0);
    }
    if (m_use_pressure_gradient_sensor && (m_sensor_source == "velocity")) {
        m_grad_vel = this->m_sim.repo().create_scratch_field(
            AMREX_SPACEDIM * AMREX_SPACEDIM, 0);
    }
    // A restart reads the gate from the checkpoint file before this point
    if ((m_separation_gate != nullptr) &&
        !this->m_sim.io_manager().is_restart()) {
        m_separation_gate->setVal(0.0_rt);
    }
}

template <typename Transport>
void KLAxellSeparation<Transport>::post_regrid_actions()
{
    KLAxell<Transport>::post_regrid_actions();
    m_geometry = this->m_sim.repo().create_scratch_field(
        klaxell_separation::geom_ncomp, 1);
    // The curvature Rt is read by update_alphaeff before the first viscosity
    // update
    for (int lev = 0; lev < this->m_sim.repo().num_active_levels(); ++lev) {
        (*m_geometry)(lev).setVal(0.0_rt);
    }
    if (m_use_pressure_gradient_sensor) {
        m_pressure_gradient_sensor->setVal(0.0_rt);
    }
    if (m_use_realizable_cmu) {
        m_strain = this->m_sim.repo().create_scratch_field(1, 0);
    }
    if (m_use_pressure_gradient_sensor && (m_sensor_source == "velocity")) {
        m_grad_vel = this->m_sim.repo().create_scratch_field(
            AMREX_SPACEDIM * AMREX_SPACEDIM, 0);
    }
}

template <typename Transport>
void KLAxellSeparation<Transport>::update_turbulent_viscosity(
    const FieldState fstate, const DiffusionType /*unused*/)
{
    BL_PROFILE(
        "kynema-sgf::" + this->identifier() + "::update_turbulent_viscosity");

    fvm::gradient(*this->m_gradT, this->m_temperature.state(fstate));

    const auto& vel = this->m_vel.state(fstate);
    fvm::strainrate(this->m_shear_prod, vel);
    const bool use_velocity_sensor =
        m_use_pressure_gradient_sensor && (m_sensor_source == "velocity");
    if (use_velocity_sensor) {
        fvm::gradient(*m_grad_vel, vel);
    }

    const auto beta = (this->m_transport).beta();
    auto& mu_turb = this->mu_turb();
    const int nlevels = mu_turb.repo().num_active_levels();
    if (m_use_realizable_cmu) {
        // The closure turns the strain rate into the shear production; keep
        // the strain rate for the limiter
        for (int lev = 0; lev < nlevels; ++lev) {
            amrex::MultiFab::Copy(
                (*m_strain)(lev), (this->m_shear_prod)(lev), 0, 0, 1, 0);
        }
    }
    const bool has_terrain =
        this->m_sim.repo().int_field_exists("terrain_blank");
    for (int lev = 0; lev < nlevels; ++lev) {
        if (has_terrain) {
            terrain_drag_geometry(lev);
        } else {
            flat_geometry(lev);
        }
        closure(lev, fstate, *beta);
        if (use_velocity_sensor) {
            velocity_sensor(lev, fstate);
        } else if (m_use_pressure_gradient_sensor) {
            pressure_gradient_sensor(lev, fstate);
        }
        if (m_use_curvature_correction && (m_curvature_model == "richardson")) {
            curvature_richardson(lev, fstate);
        } else if (m_use_curvature_correction) {
            curvature_rotation_function(lev, fstate);
        }
        if (m_use_realizable_cmu && (m_separation_gate != nullptr)) {
            realizable_cmu_relaxed(lev);
        } else if (m_use_realizable_cmu) {
            realizable_cmu(lev);
        }
    }
    amrex::Gpu::streamSynchronize();

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

template <typename Transport>
void KLAxellSeparation<Transport>::flat_geometry(const int lev)
{
    const auto& geom = this->m_sim.repo().mesh().Geom(lev);
    const auto& problo = geom.ProbLoArray();
    const amrex::Real dz = geom.CellSize()[2];
    const auto& geom_arrs = (*m_geometry)(lev).arrays();

    amrex::ParallelFor(
        (*m_geometry)(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            geom_arrs[nbx](i, j, k, klaxell_separation::geom_height) =
                problo[2] + ((k + 0.5_rt) * dz);
            geom_arrs[nbx](i, j, k, klaxell_separation::geom_fluid_weight) =
                1.0_rt;
        });
}

template <typename Transport>
void KLAxellSeparation<Transport>::terrain_drag_geometry(const int lev)
{
    const auto& geom = this->m_sim.repo().mesh().Geom(lev);
    const auto& problo = geom.ProbLoArray();
    const amrex::Real dz = geom.CellSize()[2];
    const auto& geom_arrs = (*m_geometry)(lev).arrays();
    const auto& ht_arrs =
        this->m_sim.repo().get_field("terrain_height")(lev).const_arrays();
    const auto& blank_arrs =
        this->m_sim.repo().get_int_field("terrain_blank")(lev).const_arrays();

    amrex::ParallelFor(
        (*m_geometry)(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            geom_arrs[nbx](i, j, k, klaxell_separation::geom_height) =
                amrex::max<amrex::Real>(
                    problo[2] + ((k + 0.5_rt) * dz) - ht_arrs[nbx](i, j, k),
                    0.5_rt * dz);
            geom_arrs[nbx](i, j, k, klaxell_separation::geom_fluid_weight) =
                1.0_rt - static_cast<amrex::Real>(blank_arrs[nbx](i, j, k));
        });
}

template <typename Transport>
void KLAxellSeparation<Transport>::closure(
    const int lev, const FieldState fstate, const ScratchField& beta)
{
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        this->m_gravity[0], this->m_gravity[1], this->m_gravity[2]};
    const amrex::Real Cmu = this->m_Cmu;
    const amrex::Real Cb_stable = this->m_Cb_stable;
    const amrex::Real Cb_unstable = this->m_Cb_unstable;
    const amrex::Real Rtc = -1.0_rt;
    const amrex::Real Rtmin = -3.0_rt;
    const amrex::Real lambda = 30.0_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real surf_flux = this->m_surf_flux;
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real lengthscale_switch = this->m_meso_sponge_start;

    auto& mu_turb = this->mu_turb();
    const auto& den = this->m_rho.state(fstate);
    const auto& mu_arrs = mu_turb(lev).arrays();
    const auto& rho_arrs = den(lev).const_arrays();
    const auto& gradT_arrs = (*this->m_gradT)(lev).const_arrays();
    const auto& tlscale_arrs = (this->m_turb_lscale)(lev).arrays();
    const auto& tke_arrs = (*this->m_tke)(lev).const_arrays();
    const auto& buoy_prod_arrs = (this->m_buoy_prod)(lev).arrays();
    const auto& shear_prod_arrs = (this->m_shear_prod)(lev).arrays();
    const auto& beta_arrs = beta(lev).const_arrays();
    const auto& geom_arrs = (*m_geometry)(lev).arrays();

    // Same operations, in the same order, as the KLAxell kernels so that the
    // model reproduces KLAxell when no treatment is enabled
    amrex::ParallelFor(
        mu_turb(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real z =
                geom_arrs[nbx](i, j, k, klaxell_separation::geom_height);
            const amrex::Real fluid_weight =
                geom_arrs[nbx](i, j, k, klaxell_separation::geom_fluid_weight);
            const amrex::Real stratification =
                -((gradT_arrs[nbx](i, j, k, 0) * gravity[0]) +
                  (gradT_arrs[nbx](i, j, k, 1) * gravity[1]) +
                  (gradT_arrs[nbx](i, j, k, 2) * gravity[2])) *
                beta_arrs[nbx](i, j, k);
            const amrex::Real lscale_s =
                (lambda * kappa * z) / (lambda + (kappa * z));
            const amrex::Real lscale_b =
                Cb_stable * std::sqrt(
                                tke_arrs[nbx](i, j, k) /
                                amrex::max<amrex::Real>(stratification, tiny));
            const amrex::Real epsilon =
                utils::powi(Cmu, 3) * std::pow(tke_arrs[nbx](i, j, k), 1.5_rt) /
                (tlscale_arrs[nbx](i, j, k) + tiny);
            amrex::Real Rt = utils::powi(tke_arrs[nbx](i, j, k) / epsilon, 2) *
                             stratification;
            Rt = (Rt > Rtc) ? Rt
                            : amrex::max<amrex::Real>(
                                  Rt, Rt - (utils::powi(Rt - Rtc, 2) /
                                            (Rt + Rtmin - (2.0_rt * Rtc))));
            tlscale_arrs[nbx](i, j, k) =
                (stratification > 0)
                    ? std::sqrt(
                          utils::powi(lscale_s * lscale_b, 2) /
                          (utils::powi(lscale_s, 2) + utils::powi(lscale_b, 2)))
                    : lscale_s *
                          std::sqrt(
                              1.0_rt - (utils::powi(Cmu, 6) *
                                        utils::powi(Cb_unstable, -2) * Rt));
            tlscale_arrs[nbx](i, j, k) =
                (stratification > 0)
                    ? amrex::min<amrex::Real>(
                          tlscale_arrs[nbx](i, j, k),
                          std::sqrt(
                              Cmu * tke_arrs[nbx](i, j, k) / stratification))
                    : tlscale_arrs[nbx](i, j, k);
            tlscale_arrs[nbx](i, j, k) =
                (std::abs(surf_flux) < 1.0e-5_rt && z <= lengthscale_switch)
                    ? lscale_s
                    : tlscale_arrs[nbx](i, j, k);
            Rt = (std::abs(surf_flux) < 1.0e-5_rt && z <= lengthscale_switch)
                     ? 0.0_rt
                     : Rt;
            geom_arrs[nbx](i, j, k, klaxell_separation::geom_rt) = Rt;
            const amrex::Real Cmu_Rt =
                (Cmu + (0.108_rt * Rt)) /
                (1.0_rt + (0.308_rt * Rt) + (0.00837_rt * utils::powi(Rt, 2)));
            mu_arrs[nbx](i, j, k) =
                rho_arrs[nbx](i, j, k) * Cmu_Rt * tlscale_arrs[nbx](i, j, k) *
                std::sqrt(tke_arrs[nbx](i, j, k)) * fluid_weight;
            const amrex::Real Cmu_prime_Rt = Cmu / (1.0_rt + (0.277_rt * Rt));
            const amrex::Real muPrime = rho_arrs[nbx](i, j, k) * Cmu_prime_Rt *
                                        tlscale_arrs[nbx](i, j, k) *
                                        std::sqrt(tke_arrs[nbx](i, j, k)) *
                                        fluid_weight;
            buoy_prod_arrs[nbx](i, j, k) = -muPrime * stratification;
            shear_prod_arrs[nbx](i, j, k) *=
                shear_prod_arrs[nbx](i, j, k) * mu_arrs[nbx](i, j, k);
        });
}

template <typename Transport>
void KLAxellSeparation<Transport>::pressure_gradient_sensor(
    const int lev, const FieldState fstate)
{
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real c_u = m_sensor_velocity_weight;

    const auto& vel_arrs = this->m_vel.state(fstate)(lev).const_arrays();
    const auto& rho_arrs = this->m_rho.state(fstate)(lev).const_arrays();
    const auto& tke_arrs = (*this->m_tke)(lev).const_arrays();
    const auto& tlscale_arrs = (this->m_turb_lscale)(lev).const_arrays();
    const auto& gp_arrs =
        this->m_sim.repo().get_field("gp")(lev).const_arrays();
    const auto& geom_arrs = (*m_geometry)(lev).const_arrays();
    const auto& sensor_arrs = (*m_pressure_gradient_sensor)(lev).arrays();

    // (u / |u|) . grad(p) divided by rho (k + c_u |u|^2) / L. The divisors are
    // floored so that the sensor stays finite in still air, and the fluid
    // weight zeroes it inside the terrain.
    amrex::ParallelFor(
        (*m_pressure_gradient_sensor)(lev),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const auto& vel = vel_arrs[nbx];
            const auto& gp = gp_arrs[nbx];
            const amrex::Real umag_sqr = (vel(i, j, k, 0) * vel(i, j, k, 0)) +
                                         (vel(i, j, k, 1) * vel(i, j, k, 1)) +
                                         (vel(i, j, k, 2) * vel(i, j, k, 2));
            const amrex::Real scale =
                geom_arrs[nbx](i, j, k, klaxell_separation::geom_fluid_weight) *
                tlscale_arrs[nbx](i, j, k) /
                (amrex::max<amrex::Real>(std::sqrt(umag_sqr), tiny) *
                 amrex::max<amrex::Real>(
                     rho_arrs[nbx](i, j, k) *
                         (tke_arrs[nbx](i, j, k) + (c_u * umag_sqr)),
                     tiny));
            sensor_arrs[nbx](i, j, k) = ((vel(i, j, k, 0) * gp(i, j, k, 0)) +
                                         (vel(i, j, k, 1) * gp(i, j, k, 1)) +
                                         (vel(i, j, k, 2) * gp(i, j, k, 2))) *
                                        scale;
        });
}

template <typename Transport>
void KLAxellSeparation<Transport>::velocity_sensor(
    const int lev, const FieldState fstate)
{
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real c_u = m_sensor_velocity_weight;

    // Ghost cells of the fluid weight: copied from the neighboring boxes and
    // across periodic boundaries, 1 (fluid) on the other domain boundaries
    auto& geom_mf = (*m_geometry)(lev);
    geom_mf.setBndry(1.0_rt);
    geom_mf.FillBoundary(this->m_sim.repo().mesh().Geom(lev).periodicity());

    const auto& vel_arrs = this->m_vel.state(fstate)(lev).const_arrays();
    const auto& gradvel_arrs = (*m_grad_vel)(lev).const_arrays();
    const auto& tke_arrs = (*this->m_tke)(lev).const_arrays();
    const auto& tlscale_arrs = (this->m_turb_lscale)(lev).const_arrays();
    const auto& geom_arrs = geom_mf.const_arrays();
    const auto& sensor_arrs = (*m_pressure_gradient_sensor)(lev).arrays();

    // -(u . grad)(|u|^2 / 2) / |u| = -u_m u_n d(u_n)/d(x_m) / |u|, divided by
    // (k + c_u |u|^2) / L with the divisors floored as in the pressure sensor.
    // The gradient field stores d(u_n)/d(x_m) in component n * 3 + m. The
    // smallest fluid weight of the face neighbors zeroes the sensor next to
    // the terrain.
    amrex::ParallelFor(
        (*m_pressure_gradient_sensor)(lev),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const auto& vel = vel_arrs[nbx];
            const auto& gradvel = gradvel_arrs[nbx];
            const auto& geo = geom_arrs[nbx];
            amrex::Real advection = 0.0_rt;
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                for (int m = 0; m < AMREX_SPACEDIM; ++m) {
                    advection += vel(i, j, k, m) * vel(i, j, k, n) *
                                 gradvel(i, j, k, (n * AMREX_SPACEDIM) + m);
                }
            }
            const amrex::Real umag_sqr = (vel(i, j, k, 0) * vel(i, j, k, 0)) +
                                         (vel(i, j, k, 1) * vel(i, j, k, 1)) +
                                         (vel(i, j, k, 2) * vel(i, j, k, 2));
            const int w = klaxell_separation::geom_fluid_weight;
            const amrex::Real neighbor_weight = amrex::min(
                geo(i - 1, j, k, w), geo(i + 1, j, k, w), geo(i, j - 1, k, w),
                geo(i, j + 1, k, w), geo(i, j, k - 1, w), geo(i, j, k + 1, w));
            const amrex::Real scale =
                geo(i, j, k, w) * neighbor_weight * tlscale_arrs[nbx](i, j, k) /
                (amrex::max<amrex::Real>(std::sqrt(umag_sqr), tiny) *
                 amrex::max<amrex::Real>(
                     tke_arrs[nbx](i, j, k) + (c_u * umag_sqr), tiny));
            sensor_arrs[nbx](i, j, k) = -advection * scale;
        });
}

template <typename Transport>
void KLAxellSeparation<Transport>::realizable_cmu(const int lev)
{
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real Cmu = this->m_Cmu;
    const amrex::Real threshold = m_sensor_threshold;
    const amrex::Real strength = m_realizable_cmu_strength;

    const auto& mu_arrs = this->mu_turb()(lev).arrays();
    const auto& buoy_prod_arrs = (this->m_buoy_prod)(lev).arrays();
    const auto& shear_prod_arrs = (this->m_shear_prod)(lev).arrays();
    const auto& tke_arrs = (*this->m_tke)(lev).const_arrays();
    const auto& tlscale_arrs = (this->m_turb_lscale)(lev).const_arrays();
    const auto& strain_arrs = (*m_strain)(lev).const_arrays();
    const auto& sensor_arrs = (*m_pressure_gradient_sensor)(lev).const_arrays();

    // Divide by 1 + c_s g max(0, Sigma / Cmu - 1) with the sensor gate g; the
    // factor is exactly 1 where g = 0, so those cells keep the KLAxell values
    amrex::ParallelFor(
        this->mu_turb()(lev),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            // Linear ramp from 0 at the threshold to 1 at twice the threshold,
            // exactly 0 below it
            const amrex::Real gate = amrex::min<amrex::Real>(
                amrex::max<amrex::Real>(
                    (sensor_arrs[nbx](i, j, k) - threshold) / threshold,
                    0.0_rt),
                1.0_rt);
            const amrex::Real sigma =
                tlscale_arrs[nbx](i, j, k) * strain_arrs[nbx](i, j, k) /
                amrex::max<amrex::Real>(
                    std::sqrt(tke_arrs[nbx](i, j, k)), tiny);
            const amrex::Real factor =
                1.0_rt / (1.0_rt + (strength * gate *
                                    amrex::max<amrex::Real>(
                                        (sigma / Cmu) - 1.0_rt, 0.0_rt)));
            mu_arrs[nbx](i, j, k) *= factor;
            buoy_prod_arrs[nbx](i, j, k) *= factor;
            shear_prod_arrs[nbx](i, j, k) *= factor;
        });
}

template <typename Transport>
void KLAxellSeparation<Transport>::realizable_cmu_relaxed(const int lev)
{
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real Cmu = this->m_Cmu;
    const amrex::Real strength = m_realizable_cmu_strength;

    const auto& mu_arrs = this->mu_turb()(lev).arrays();
    const auto& buoy_prod_arrs = (this->m_buoy_prod)(lev).arrays();
    const auto& shear_prod_arrs = (this->m_shear_prod)(lev).arrays();
    const auto& tke_arrs = (*this->m_tke)(lev).const_arrays();
    const auto& tlscale_arrs = (this->m_turb_lscale)(lev).const_arrays();
    const auto& strain_arrs = (*m_strain)(lev).const_arrays();
    const auto& gate_arrs = (*m_separation_gate)(lev).const_arrays();

    // Same factor as realizable_cmu with the stored gate
    amrex::ParallelFor(
        this->mu_turb()(lev),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real sigma =
                tlscale_arrs[nbx](i, j, k) * strain_arrs[nbx](i, j, k) /
                amrex::max<amrex::Real>(
                    std::sqrt(tke_arrs[nbx](i, j, k)), tiny);
            const amrex::Real factor =
                1.0_rt / (1.0_rt + (strength * gate_arrs[nbx](i, j, k) *
                                    amrex::max<amrex::Real>(
                                        (sigma / Cmu) - 1.0_rt, 0.0_rt)));
            mu_arrs[nbx](i, j, k) *= factor;
            buoy_prod_arrs[nbx](i, j, k) *= factor;
            shear_prod_arrs[nbx](i, j, k) *= factor;
        });
}

template <typename Transport>
void KLAxellSeparation<Transport>::post_advance_work()
{
    KLAxell<Transport>::post_advance_work();
    if (m_separation_gate == nullptr) {
        return;
    }
    const int nlevels = this->m_sim.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        relax_gate(lev);
    }
    amrex::Gpu::streamSynchronize();
}

template <typename Transport>
void KLAxellSeparation<Transport>::relax_gate(const int lev)
{
    const amrex::Real threshold = m_sensor_threshold;
    const amrex::Real decay =
        std::exp(-this->m_sim.time().delta_t() / m_gate_relaxation_time);

    const auto& sensor_arrs = (*m_pressure_gradient_sensor)(lev).const_arrays();
    const auto& gate_arrs = (*m_separation_gate)(lev).arrays();

    // The sensor is the one of the last viscosity update of the step
    amrex::ParallelFor(
        (*m_separation_gate)(lev),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real target = amrex::min<amrex::Real>(
                amrex::max<amrex::Real>(
                    (sensor_arrs[nbx](i, j, k) - threshold) / threshold,
                    0.0_rt),
                1.0_rt);
            gate_arrs[nbx](i, j, k) =
                target + ((gate_arrs[nbx](i, j, k) - target) * decay);
        });
}

template <typename Transport>
void KLAxellSeparation<Transport>::curvature_richardson(
    const int lev, const FieldState fstate)
{
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real Cmu = this->m_Cmu;
    const amrex::Real Cmu6 = utils::powi(Cmu, 6);
    const amrex::Real coeff = m_curvature_coefficient;
    const amrex::Real tolerance = klaxell_separation::curvature_tolerance;
    const amrex::Real Rtc = -1.0_rt;
    const amrex::Real Rtmin = -3.0_rt;
    const auto& geom = this->m_sim.repo().mesh().Geom(lev);
    const auto idx = geom.InvCellSizeArray();

    // Ghost cells of the fluid weight, as in velocity_sensor
    auto& geom_mf = (*m_geometry)(lev);
    geom_mf.setBndry(1.0_rt);
    geom_mf.FillBoundary(geom.periodicity());

    const auto& vel_arrs = this->m_vel.state(fstate)(lev).const_arrays();
    const auto& tke_arrs = (*this->m_tke)(lev).const_arrays();
    const auto& tlscale_arrs = (this->m_turb_lscale)(lev).const_arrays();
    const auto& mu_arrs = this->mu_turb()(lev).arrays();
    const auto& buoy_prod_arrs = (this->m_buoy_prod)(lev).arrays();
    const auto& shear_prod_arrs = (this->m_shear_prod)(lev).arrays();
    const auto& geom_arrs = geom_mf.arrays();

    // N2_c = 2 (|a_n|^2 - a_n . grad(|u|^2 / 2)) / |u|^2 with a = (u . grad) u
    // and a_n = a - (a . u / |u|^2) u, from central differences. Below the
    // tolerance (relative to the squared strain rate) the cell is unchanged.
    amrex::ParallelFor(
        this->mu_turb()(lev),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const auto& vel = vel_arrs[nbx];
            const auto& geo = geom_arrs[nbx];
            // grad[n][m] = d(u_n) / d(x_m)
            amrex::GpuArray<
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>, AMREX_SPACEDIM>
                grad{};
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                grad[n][0] = 0.5_rt *
                             (vel(i + 1, j, k, n) - vel(i - 1, j, k, n)) *
                             idx[0];
                grad[n][1] = 0.5_rt *
                             (vel(i, j + 1, k, n) - vel(i, j - 1, k, n)) *
                             idx[1];
                grad[n][2] = 0.5_rt *
                             (vel(i, j, k + 1, n) - vel(i, j, k - 1, n)) *
                             idx[2];
            }
            amrex::Real umag_sqr = 0.0_rt;
            amrex::Real strain_sqr = 0.0_rt;
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> accel{};
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> grad_q{};
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                umag_sqr += vel(i, j, k, n) * vel(i, j, k, n);
                accel[n] = 0.0_rt;
                grad_q[n] = 0.0_rt;
            }
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                for (int m = 0; m < AMREX_SPACEDIM; ++m) {
                    accel[n] += vel(i, j, k, m) * grad[n][m];
                    grad_q[m] += vel(i, j, k, n) * grad[n][m];
                    const amrex::Real s = 0.5_rt * (grad[n][m] + grad[m][n]);
                    strain_sqr += 2.0_rt * s * s;
                }
            }
            const amrex::Real inv_umag_sqr =
                1.0_rt / amrex::max<amrex::Real>(umag_sqr, tiny);
            amrex::Real a_dot_u = 0.0_rt;
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                a_dot_u += accel[n] * vel(i, j, k, n);
            }
            amrex::Real a_n_sqr = 0.0_rt;
            amrex::Real a_n_dot_grad_q = 0.0_rt;
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                const amrex::Real a_n =
                    accel[n] - (a_dot_u * inv_umag_sqr * vel(i, j, k, n));
                a_n_sqr += a_n * a_n;
                a_n_dot_grad_q += a_n * grad_q[n];
            }
            const amrex::Real n2_curv =
                2.0_rt * (a_n_sqr - a_n_dot_grad_q) * inv_umag_sqr;

            const int w = klaxell_separation::geom_fluid_weight;
            const amrex::Real weight =
                geo(i, j, k, w) * amrex::min(
                                      geo(i - 1, j, k, w), geo(i + 1, j, k, w),
                                      geo(i, j - 1, k, w), geo(i, j + 1, k, w),
                                      geo(i, j, k - 1, w), geo(i, j, k + 1, w));
            const bool active =
                (weight > 0.0_rt) &&
                (std::abs(n2_curv) >
                 tolerance * amrex::max<amrex::Real>(strain_sqr, tiny));
            const amrex::Real lscale = tlscale_arrs[nbx](i, j, k);
            const amrex::Real rt_curv =
                active ? coeff * weight * lscale * lscale * n2_curv /
                             (Cmu6 * amrex::max<amrex::Real>(
                                         tke_arrs[nbx](i, j, k), tiny))
                       : 0.0_rt;

            const amrex::Real rt = geo(i, j, k, klaxell_separation::geom_rt);
            amrex::Real rt_tot = rt + rt_curv;
            rt_tot =
                (rt_tot > Rtc)
                    ? rt_tot
                    : amrex::max<amrex::Real>(
                          rt_tot, rt_tot - (utils::powi(rt_tot - Rtc, 2) /
                                            (rt_tot + Rtmin - (2.0_rt * Rtc))));
            const amrex::Real cmu_rt =
                (Cmu + (0.108_rt * rt)) /
                (1.0_rt + (0.308_rt * rt) + (0.00837_rt * rt * rt));
            const amrex::Real cmu_rt_tot =
                (Cmu + (0.108_rt * rt_tot)) /
                (1.0_rt + (0.308_rt * rt_tot) + (0.00837_rt * rt_tot * rt_tot));
            const amrex::Real mu_factor = active ? cmu_rt_tot / cmu_rt : 1.0_rt;
            const amrex::Real prime_factor =
                active ? (1.0_rt + (0.277_rt * rt)) /
                             (1.0_rt + (0.277_rt * rt_tot))
                       : 1.0_rt;
            mu_arrs[nbx](i, j, k) *= mu_factor;
            shear_prod_arrs[nbx](i, j, k) *= mu_factor;
            buoy_prod_arrs[nbx](i, j, k) *= prime_factor;
            geo(i, j, k, klaxell_separation::geom_rt_curvature) = rt_curv;
        });
}

template <typename Transport>
void KLAxellSeparation<Transport>::curvature_rotation_function(
    const int lev, const FieldState fstate)
{
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real cr1 = m_curvature_cr1;
    const amrex::Real cr2 = m_curvature_cr2;
    const amrex::Real cr3 = m_curvature_cr3;
    const amrex::Real tolerance =
        klaxell_separation::rotation_function_tolerance;
    const auto& geom = this->m_sim.repo().mesh().Geom(lev);
    const auto idx = geom.InvCellSizeArray();

    // Ghost cells of the fluid weight, as in velocity_sensor
    auto& geom_mf = (*m_geometry)(lev);
    geom_mf.setBndry(1.0_rt);
    geom_mf.FillBoundary(geom.periodicity());

    const auto& vel_arrs = this->m_vel.state(fstate)(lev).const_arrays();
    const auto& mu_arrs = this->mu_turb()(lev).arrays();
    const auto& buoy_prod_arrs = (this->m_buoy_prod)(lev).arrays();
    const auto& shear_prod_arrs = (this->m_shear_prod)(lev).arrays();
    const auto& geom_arrs = geom_mf.const_arrays();

    // First and second velocity derivatives by central differences; the
    // cell is unchanged where the flow has no strain or rotation or where f
    // stays within the tolerance of 1
    amrex::ParallelFor(
        this->mu_turb()(lev),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const auto& vel = vel_arrs[nbx];
            const auto& geo = geom_arrs[nbx];
            // grad[n][a] = d(u_n)/d(x_a), hess[n][a][b] = d2(u_n)/d(x_a)d(x_b)
            amrex::GpuArray<
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>, AMREX_SPACEDIM>
                grad{};
            amrex::GpuArray<
                amrex::GpuArray<
                    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>,
                    AMREX_SPACEDIM>,
                AMREX_SPACEDIM>
                hess{};
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                for (int a = 0; a < AMREX_SPACEDIM; ++a) {
                    const int ai = static_cast<int>(a == 0);
                    const int aj = static_cast<int>(a == 1);
                    const int ak = static_cast<int>(a == 2);
                    grad[n][a] = 0.5_rt *
                                 (vel(i + ai, j + aj, k + ak, n) -
                                  vel(i - ai, j - aj, k - ak, n)) *
                                 idx[a];
                    for (int b = 0; b < AMREX_SPACEDIM; ++b) {
                        const int bi = static_cast<int>(b == 0);
                        const int bj = static_cast<int>(b == 1);
                        const int bk = static_cast<int>(b == 2);
                        hess[n][a][b] = (a == b)
                                            ? (vel(i + ai, j + aj, k + ak, n) -
                                               (2.0_rt * vel(i, j, k, n)) +
                                               vel(i - ai, j - aj, k - ak, n)) *
                                                  idx[a] * idx[a]
                                            : 0.25_rt *
                                                  (vel(i + ai + bi, j + aj + bj,
                                                       k + ak + bk, n) -
                                                   vel(i + ai - bi, j + aj - bj,
                                                       k + ak - bk, n) -
                                                   vel(i - ai + bi, j - aj + bj,
                                                       k - ak + bk, n) +
                                                   vel(i - ai - bi, j - aj - bj,
                                                       k - ak - bk, n)) *
                                                  idx[a] * idx[b];
                    }
                }
            }
            // Strain and rotation tensors, their magnitudes and the material
            // derivative of the strain, u_c d(S_pq)/d(x_c)
            amrex::GpuArray<
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>, AMREX_SPACEDIM>
                strain{};
            amrex::GpuArray<
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>, AMREX_SPACEDIM>
                rotation{};
            amrex::GpuArray<
                amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>, AMREX_SPACEDIM>
                dstrain{};
            amrex::Real strain_sqr = 0.0_rt;
            amrex::Real rotation_sqr = 0.0_rt;
            for (int p = 0; p < AMREX_SPACEDIM; ++p) {
                for (int q = 0; q < AMREX_SPACEDIM; ++q) {
                    strain[p][q] = 0.5_rt * (grad[p][q] + grad[q][p]);
                    rotation[p][q] = 0.5_rt * (grad[p][q] - grad[q][p]);
                    strain_sqr += 2.0_rt * strain[p][q] * strain[p][q];
                    rotation_sqr += 2.0_rt * rotation[p][q] * rotation[p][q];
                    dstrain[p][q] = 0.0_rt;
                    for (int c = 0; c < AMREX_SPACEDIM; ++c) {
                        dstrain[p][q] += 0.5_rt * vel(i, j, k, c) *
                                         (hess[p][q][c] + hess[q][p][c]);
                    }
                }
            }
            const amrex::Real d_sqr = 0.5_rt * (strain_sqr + rotation_sqr);
            amrex::Real r_tilde_num = 0.0_rt;
            for (int p = 0; p < AMREX_SPACEDIM; ++p) {
                for (int q = 0; q < AMREX_SPACEDIM; ++q) {
                    for (int c = 0; c < AMREX_SPACEDIM; ++c) {
                        r_tilde_num += 2.0_rt * rotation[p][c] * strain[q][c] *
                                       dstrain[p][q];
                    }
                }
            }
            const amrex::Real r_star =
                std::sqrt(strain_sqr) /
                amrex::max<amrex::Real>(std::sqrt(rotation_sqr), tiny);
            const amrex::Real r_tilde =
                r_tilde_num / amrex::max<amrex::Real>(d_sqr * d_sqr, tiny);
            const amrex::Real f_rot = amrex::min<amrex::Real>(
                amrex::max<amrex::Real>(
                    ((1.0_rt + cr1) * (2.0_rt * r_star / (1.0_rt + r_star)) *
                     (1.0_rt - (cr3 * std::atan(cr2 * r_tilde)))) -
                        cr1,
                    0.0_rt),
                1.25_rt);

            const int w = klaxell_separation::geom_fluid_weight;
            const amrex::Real weight =
                geo(i, j, k, w) * amrex::min(
                                      geo(i - 1, j, k, w), geo(i + 1, j, k, w),
                                      geo(i, j - 1, k, w), geo(i, j + 1, k, w),
                                      geo(i, j, k - 1, w), geo(i, j, k + 1, w));
            // Deviations of f from 1 below the tolerance are ignored; above it
            // they ramp in linearly to the full deviation at twice the
            // tolerance
            const amrex::Real deviation = f_rot - 1.0_rt;
            const amrex::Real ramp = amrex::min<amrex::Real>(
                amrex::max<amrex::Real>(
                    (std::abs(deviation) - tolerance) / tolerance, 0.0_rt),
                1.0_rt);
            const bool active =
                (weight > 0.0_rt) && (d_sqr > tiny) && (ramp > 0.0_rt);
            const amrex::Real factor =
                active ? 1.0_rt + (weight * ramp * deviation) : 1.0_rt;
            mu_arrs[nbx](i, j, k) *= factor;
            shear_prod_arrs[nbx](i, j, k) *= factor;
            buoy_prod_arrs[nbx](i, j, k) *= factor;
        });
}

template <typename Transport>
void KLAxellSeparation<Transport>::update_alphaeff(Field& alphaeff)
{
    if (m_use_curvature_correction && (m_curvature_model == "richardson")) {
        update_alphaeff_curvature(alphaeff);
    } else {
        KLAxell<Transport>::update_alphaeff(alphaeff);
    }
}

template <typename Transport>
void KLAxellSeparation<Transport>::update_alphaeff_curvature(Field& alphaeff)
{
    BL_PROFILE(
        "kynema-sgf::" + this->identifier() + "::update_alphaeff_curvature");
    auto lam_alpha = (this->m_transport).alpha();
    auto& mu_turb = this->mu_turb();
    auto& repo = mu_turb.repo();

    fvm::gradient(*this->m_gradT, this->m_temperature);
    auto& gradT = *this->m_gradT;
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        this->m_gravity[0], this->m_gravity[1], this->m_gravity[2]};
    const auto beta = (this->m_transport).beta();
    const amrex::Real Cmu = this->m_Cmu;
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& muturb_arrs = mu_turb(lev).arrays();
        const auto& alphaeff_arrs = alphaeff(lev).arrays();
        const auto& lam_diff_arrs = (*lam_alpha)(lev).arrays();
        const auto& tke_arrs = (*this->m_tke)(lev).arrays();
        const auto& gradT_arrs = gradT(lev).const_arrays();
        const auto& tlscale_arrs = (this->m_turb_lscale)(lev).arrays();
        const auto& beta_arrs = (*beta)(lev).const_arrays();
        const auto& geom_arrs = (*m_geometry)(lev).const_arrays();
        const amrex::Real Rtc = -1.0_rt;
        const amrex::Real Rtmin = -3.0_rt;
        // Same operations as KLAxell::update_alphaeff, with the stored
        // curvature contribution added to Rt before it is limited
        amrex::ParallelFor(
            mu_turb(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                amrex::Real stratification =
                    -((gradT_arrs[nbx](i, j, k, 0) * gravity[0]) +
                      (gradT_arrs[nbx](i, j, k, 1) * gravity[1]) +
                      (gradT_arrs[nbx](i, j, k, 2) * gravity[2])) *
                    beta_arrs[nbx](i, j, k);
                amrex::Real epsilon = utils::powi(Cmu, 3) *
                                      std::pow(tke_arrs[nbx](i, j, k), 1.5_rt) /
                                      tlscale_arrs[nbx](i, j, k);
                amrex::Real Rt =
                    utils::powi(tke_arrs[nbx](i, j, k) / epsilon, 2) *
                    stratification;
                Rt += geom_arrs[nbx](
                    i, j, k, klaxell_separation::geom_rt_curvature);
                Rt = (Rt > Rtc) ? Rt
                                : amrex::max<amrex::Real>(
                                      Rt, Rt - (utils::powi(Rt - Rtc, 2) /
                                                (Rt + Rtmin - 2.0_rt * Rtc)));
                const amrex::Real prandtlRt =
                    (1.0_rt + 0.193_rt * Rt) / (1.0_rt + 0.0302_rt * Rt);
                alphaeff_arrs[nbx](i, j, k) =
                    lam_diff_arrs[nbx](i, j, k) +
                    (muturb_arrs[nbx](i, j, k) / prandtlRt);
            });
    }
    amrex::Gpu::streamSynchronize();

    alphaeff.fillpatch(this->m_sim.time().current_time());
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(KLAxellSeparation);

} // namespace kynema_sgf
