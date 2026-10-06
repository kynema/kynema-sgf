#include "src/turbulence/RANS/KLAxell.H"
#include "src/equation_systems/PDEBase.H"
#include "src/turbulence/TurbModelDefs.H"
#include "src/fvm/gradient.H"
#include "src/fvm/strainrate.H"
#include "src/turbulence/turb_utils.H"
#include "src/equation_systems/tke/TKE.H"
#include "AMReX_ParmParse.H"
#include "src/utilities/math_ops.H"
#include "src/wind_energy/MOData.H"
#include "src/physics/TerrainDrag.H"

using namespace amrex::literals;

namespace kynema_sgf {
namespace turbulence {

namespace {
//! Derivative at the first of three points p0, p1, p2 spaced s and t along
//! the direction away from a wall, from the quadratic through them. For
//! s = t = dx it is the stencil of flat ground at its bottom wall for a
//! tangential velocity: the zlo stencil (p1 / 3 + p0 - 4 g / 3) / dx with the
//! hoextrap ghost g = (15 p0 - 10 p1 + 3 p2) / 8, i.e.
//! (-3 p0 + 4 p1 - p2) / (2 dx).
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real one_sided_derivative(
    const amrex::Real p0,
    const amrex::Real p1,
    const amrex::Real p2,
    const amrex::Real s,
    const amrex::Real t)
{
    return (-((2.0_rt * s) + t) / (s * (s + t)) * p0) +
           ((s + t) / (s * t) * p1) - (s / ((s + t) * t) * p2);
}

//! Derivative at distance a from a wall where the value is zero, from the
//! quadratic through the wall, the point p0 and the next point p1 at a + s.
//! For a = dx / 2 and s = dx it is the stencil of flat ground at its bottom
//! wall for the wall-normal velocity, whose ghost holds the zero wall value:
//! (p1 / 3 + p0) / dx.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real zero_wall_derivative(
    const amrex::Real p0,
    const amrex::Real p1,
    const amrex::Real a,
    const amrex::Real s)
{
    return ((s - a) / (a * s) * p0) + (a / ((a + s) * s) * p1);
}

//! Strain-rate magnitude sqrt(2 S_ij S_ij) of cell (i, j, k) as
//! fvm::strainrate computes it, except toward terrain. side[2 d] (low) and
//! side[2 d + 1] (high) are 0 for a fluid neighbor, 1 for a blanked cell and
//! 2 for a non-periodic domain face. Toward a blanked cell the derivative
//! takes the stencil flat ground uses at its bottom wall, the wall being the
//! face between the two cells: one-sided over the fluid cells for a
//! tangential component, through the zero wall value for the normal one.
//! Toward a domain face it keeps the boundary stencil of fvm::strainrate;
//! between two blanked cells the derivative is zero.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real wall_strain_rate(
    amrex::Array4<amrex::Real const> const& vel,
    const int i,
    const int j,
    const int k,
    const amrex::GpuArray<int, 2 * AMREX_SPACEDIM>& side,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx)
{
    // g[c][d] = d u_c / d x_d
    amrex::Real g[AMREX_SPACEDIM][AMREX_SPACEDIM];
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const int di = (d == 0) ? 1 : 0;
        const int dj = (d == 1) ? 1 : 0;
        const int dk = (d == 2) ? 1 : 0;
        const int lo = side[2 * d];
        const int hi = side[(2 * d) + 1];
        for (int c = 0; c < AMREX_SPACEDIM; ++c) {
            const amrex::Real p0 = vel(i, j, k, c);
            const amrex::Real pp = vel(i + di, j + dj, k + dk, c);
            const amrex::Real pm = vel(i - di, j - dj, k - dk, c);
            amrex::Real grad = 0.0_rt;
            if (lo == 2) {
                grad = ((pp / 3.0_rt) + p0 - (4.0_rt / 3.0_rt * pm)) / dx[d];
            } else if (hi == 2) {
                grad = ((4.0_rt / 3.0_rt * pp) - p0 - (pm / 3.0_rt)) / dx[d];
            } else if ((lo == 1) && (hi == 1)) {
                grad = 0.0_rt;
            } else if (lo == 1) {
                grad =
                    (c == d)
                        ? zero_wall_derivative(p0, pp, 0.5_rt * dx[d], dx[d])
                        : one_sided_derivative(
                              p0, pp,
                              vel(i + (2 * di), j + (2 * dj), k + (2 * dk), c),
                              dx[d], dx[d]);
            } else if (hi == 1) {
                grad = -(
                    (c == d)
                        ? zero_wall_derivative(p0, pm, 0.5_rt * dx[d], dx[d])
                        : one_sided_derivative(
                              p0, pm,
                              vel(i - (2 * di), j - (2 * dj), k - (2 * dk), c),
                              dx[d], dx[d]));
            } else {
                grad = 0.5_rt * (pp - pm) / dx[d];
            }
            g[c][d] = grad;
        }
    }
    return std::sqrt(
        (2.0_rt * g[0][0] * g[0][0]) + (2.0_rt * g[1][1] * g[1][1]) +
        (2.0_rt * g[2][2] * g[2][2]) +
        ((g[0][1] + g[1][0]) * (g[0][1] + g[1][0])) +
        ((g[1][2] + g[2][1]) * (g[1][2] + g[2][1])) +
        ((g[2][0] + g[0][2]) * (g[2][0] + g[0][2])));
}
} // namespace

template <typename Transport>
KLAxell<Transport>::KLAxell(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_turb_lscale(sim.repo().declare_field("turb_lscale", 1))
    , m_shear_prod(sim.repo().declare_field("shear_prod", 1))
    , m_buoy_prod(sim.repo().declare_field("buoy_prod", 1))
    , m_dissip(sim.repo().declare_field("dissipation", 1))
    , m_rho(sim.repo().get_field("density"))
    , m_temperature(sim.repo().get_field("temperature"))
{
    auto& tke_eqn =
        sim.pde_manager().register_transport_pde(pde::TKE::pde_name());
    m_tke = &(tke_eqn.fields().field);
    auto& phy_mgr = this->m_sim.physics_manager();
    if (!phy_mgr.contains("ABL")) {
        amrex::Abort("KLAxell model only works with ABL physics");
    }
    {
        amrex::ParmParse pp("ABL");
        pp.get("surface_temp_flux", m_surf_flux);
        pp.query("meso_sponge_start", m_meso_sponge_start);
    }

    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }

    {
        // Near-wall options for TerrainDrag: all off with the original
        // wall treatment (default), all on with the improved one, and each
        // can be set on its own
        const bool improved = terraindrag::improved_wall_treatment();
        m_terrain_wall_stencil = improved;
        m_terrain_blanked_face_length = improved;
        m_terrain_face_stress = improved;
        m_terrain_face_heat_flux = improved;
        amrex::ParmParse pp("KLAxell");
        pp.query("terrain_wall_stencil", m_terrain_wall_stencil);
        pp.query("terrain_blanked_face_length", m_terrain_blanked_face_length);
        pp.query("terrain_face_stress", m_terrain_face_stress);
        pp.query("terrain_face_heat_flux", m_terrain_face_heat_flux);
    }
    {
        amrex::ParmParse pp_abl("ABL");
        pp_abl.query("wall_het_model", m_wall_het_model);
        pp_abl.query("monin_obukhov_length", m_monin_obukhov_length);
        pp_abl.query("kappa", m_kappa);
        pp_abl.query("mo_beta_m", m_beta_m);
        pp_abl.query("mo_gamma_m", m_gamma_m);
        amrex::ParmParse pp_drag("DragForcing");
        pp_drag.query("minimum_z0", m_min_z0);
    }

    // TKE source term to be added to PDE
    turb_utils::inject_turbulence_src_terms(
        pde::TKE::pde_name(), {"KransAxell"});
}

template <typename Transport>
void KLAxell<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Cmu", this->m_Cmu);
    pp.query("Cmu_prime", this->m_Cmu_prime);
    pp.query("Cb_stable", this->m_Cb_stable);
    pp.query("Cb_unstable", this->m_Cb_unstable);
    pp.query("prandtl", this->m_prandtl);
}

template <typename Transport>
TurbulenceModel::CoeffsDictType KLAxell<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{
        {"Cmu", this->m_Cmu},
        {"Cmu_prime", this->m_Cmu_prime},
        {"Cb_stable", this->m_Cb_stable},
        {"Cb_unstable", this->m_Cb_unstable},
        {"prandtl", this->m_prandtl}};
}

template <typename Transport>
void KLAxell<Transport>::post_init_actions()
{
    m_gradT = this->m_sim.repo().create_scratch_field(3, 0);
}

template <typename Transport>
void KLAxell<Transport>::post_regrid_actions()
{
    m_gradT = this->m_sim.repo().create_scratch_field(3, 0);
}

template <typename Transport>
void KLAxell<Transport>::update_turbulent_viscosity(
    const FieldState fstate, const DiffusionType /*unused*/)
{
    BL_PROFILE(
        "kynema-sgf::" + this->identifier() + "::update_turbulent_viscosity");

    fvm::gradient(*m_gradT, m_temperature.state(fstate));
    auto& gradT = *m_gradT;

    const auto& vel = this->m_vel.state(fstate);
    fvm::strainrate(this->m_shear_prod, vel);
    if (m_terrain_wall_stencil &&
        this->m_sim.repo().int_field_exists("terrain_blank")) {
        terrain_wall_strain_rate(vel);
    }

    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};
    const auto beta = (this->m_transport).beta();
    const amrex::Real Cmu = m_Cmu;
    const amrex::Real Cb_stable = m_Cb_stable;
    const amrex::Real Cb_unstable = m_Cb_unstable;
    auto& mu_turb = this->mu_turb();
    const auto& den = this->m_rho.state(fstate);
    const auto& repo = mu_turb.repo();
    const auto& geom_vec = repo.mesh().Geom();
    const int nlevels = repo.num_active_levels();
    const amrex::Real Rtc = -1.0_rt;
    const amrex::Real Rtmin = -3.0_rt;
    const amrex::Real lambda = 30.0_rt;
    const amrex::Real kappa = 0.41_rt;
    const amrex::Real surf_flux = m_surf_flux;
    const auto tiny = std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real lengthscale_switch = m_meso_sponge_start;
    // KLAxell.terrain_blanked_face_length: the TerrainDrag kernel below reads
    // the height of the top face of the blanked column in place of the
    // terrain height
    std::unique_ptr<ScratchField> face_height;
    if (m_terrain_blanked_face_length &&
        this->m_sim.repo().int_field_exists("terrain_blank")) {
        face_height = repo.create_scratch_field(1, 0);
        terrain_blanked_face_height(*face_height);
    }
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];
        const auto& problo = repo.mesh().Geom(lev).ProbLoArray();
        const amrex::Real dz = geom.CellSize()[2];

        const auto& mu_arrs = mu_turb(lev).arrays();
        const auto& rho_arrs = den(lev).const_arrays();
        const auto& gradT_arrs = gradT(lev).const_arrays();
        const auto& tlscale_arrs = (this->m_turb_lscale)(lev).arrays();
        const auto& tke_arrs = (*this->m_tke)(lev).arrays();
        const auto& buoy_prod_arrs = (this->m_buoy_prod)(lev).arrays();
        const auto& shear_prod_arrs = (this->m_shear_prod)(lev).arrays();
        const auto& beta_arrs = (*beta)(lev).const_arrays();

        //! Add terrain components
        const bool has_terrain =
            this->m_sim.repo().int_field_exists("terrain_blank");
        if (has_terrain) {
            const auto* m_terrain_height =
                &this->m_sim.repo().get_field("terrain_height");
            const auto* m_terrain_blank =
                &this->m_sim.repo().get_int_field("terrain_blank");
            const auto& ht_arrs = face_height
                                      ? (*face_height)(lev).const_arrays()
                                      : (*m_terrain_height)(lev).const_arrays();
            const auto& blank_arrs = (*m_terrain_blank)(lev).const_arrays();
            amrex::ParallelFor(
                mu_turb(lev),
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    amrex::Real stratification =
                        -((gradT_arrs[nbx](i, j, k, 0) * gravity[0]) +
                          (gradT_arrs[nbx](i, j, k, 1) * gravity[1]) +
                          (gradT_arrs[nbx](i, j, k, 2) * gravity[2])) *
                        beta_arrs[nbx](i, j, k);
                    const amrex::Real z = amrex::max<amrex::Real>(
                        problo[2] + ((k + 0.5_rt) * dz) - ht_arrs[nbx](i, j, k),
                        0.5_rt * dz);
                    const amrex::Real lscale_s =
                        (lambda * kappa * z) / (lambda + (kappa * z));
                    const amrex::Real lscale_b =
                        Cb_stable *
                        std::sqrt(
                            tke_arrs[nbx](i, j, k) /
                            amrex::max<amrex::Real>(stratification, tiny));
                    amrex::Real epsilon =
                        utils::powi(Cmu, 3) *
                        std::pow(tke_arrs[nbx](i, j, k), 1.5_rt) /
                        (tlscale_arrs[nbx](i, j, k) + tiny);
                    amrex::Real Rt =
                        utils::powi(tke_arrs[nbx](i, j, k) / epsilon, 2) *
                        stratification;
                    Rt = (Rt > Rtc)
                             ? Rt
                             : amrex::max<amrex::Real>(
                                   Rt, Rt - (utils::powi(Rt - Rtc, 2) /
                                             (Rt + Rtmin - (2.0_rt * Rtc))));
                    tlscale_arrs[nbx](i, j, k) =
                        (stratification > 0)
                            ? std::sqrt(
                                  utils::powi(lscale_s * lscale_b, 2) /
                                  (utils::powi(lscale_s, 2) +
                                   utils::powi(lscale_b, 2)))
                            : lscale_s *
                                  std::sqrt(
                                      1.0_rt -
                                      (utils::powi(Cmu, 6) *
                                       utils::powi(Cb_unstable, -2) * Rt));
                    tlscale_arrs[nbx](i, j, k) =
                        (stratification > 0)
                            ? amrex::min<amrex::Real>(
                                  tlscale_arrs[nbx](i, j, k),
                                  std::sqrt(
                                      Cmu * tke_arrs[nbx](i, j, k) /
                                      stratification))
                            : tlscale_arrs[nbx](i, j, k);
                    tlscale_arrs[nbx](i, j, k) =
                        (std::abs(surf_flux) < 1.0e-5_rt &&
                         z <= lengthscale_switch)
                            ? lscale_s
                            : tlscale_arrs[nbx](i, j, k);
                    Rt = (std::abs(surf_flux) < 1.0e-5_rt &&
                          z <= lengthscale_switch)
                             ? 0.0_rt
                             : Rt;
                    const amrex::Real Cmu_Rt =
                        (Cmu + (0.108_rt * Rt)) /
                        (1.0_rt + (0.308_rt * Rt) +
                         (0.00837_rt * utils::powi(Rt, 2)));
                    mu_arrs[nbx](i, j, k) = rho_arrs[nbx](i, j, k) * Cmu_Rt *
                                            tlscale_arrs[nbx](i, j, k) *
                                            std::sqrt(tke_arrs[nbx](i, j, k)) *
                                            (1.0_rt - blank_arrs[nbx](i, j, k));
                    const amrex::Real Cmu_prime_Rt =
                        Cmu / (1.0_rt + (0.277_rt * Rt));
                    const amrex::Real muPrime =
                        rho_arrs[nbx](i, j, k) * Cmu_prime_Rt *
                        tlscale_arrs[nbx](i, j, k) *
                        std::sqrt(tke_arrs[nbx](i, j, k)) *
                        (1.0_rt - blank_arrs[nbx](i, j, k));
                    buoy_prod_arrs[nbx](i, j, k) = -muPrime * stratification;
                    shear_prod_arrs[nbx](i, j, k) *=
                        shear_prod_arrs[nbx](i, j, k) * mu_arrs[nbx](i, j, k);
                });
            if (m_terrain_face_stress) {
                terrain_face_stress(lev, vel, den);
            }
        } else {
            amrex::ParallelFor(
                mu_turb(lev),
                [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    amrex::Real stratification =
                        -((gradT_arrs[nbx](i, j, k, 0) * gravity[0]) +
                          (gradT_arrs[nbx](i, j, k, 1) * gravity[1]) +
                          (gradT_arrs[nbx](i, j, k, 2) * gravity[2])) *
                        beta_arrs[nbx](i, j, k);
                    const amrex::Real z = problo[2] + ((k + 0.5_rt) * dz);
                    const amrex::Real lscale_s =
                        (lambda * kappa * z) / (lambda + (kappa * z));
                    const amrex::Real lscale_b =
                        Cb_stable *
                        std::sqrt(
                            tke_arrs[nbx](i, j, k) /
                            amrex::max<amrex::Real>(stratification, tiny));
                    amrex::Real epsilon =
                        utils::powi(Cmu, 3) *
                        std::pow(tke_arrs[nbx](i, j, k), 1.5_rt) /
                        (tlscale_arrs[nbx](i, j, k) + tiny);
                    amrex::Real Rt =
                        utils::powi(tke_arrs[nbx](i, j, k) / epsilon, 2) *
                        stratification;
                    Rt = (Rt > Rtc)
                             ? Rt
                             : amrex::max<amrex::Real>(
                                   Rt, Rt - (utils::powi(Rt - Rtc, 2) /
                                             (Rt + Rtmin - (2.0_rt * Rtc))));
                    tlscale_arrs[nbx](i, j, k) =
                        (stratification > 0)
                            ? std::sqrt(
                                  utils::powi(lscale_s * lscale_b, 2) /
                                  (utils::powi(lscale_s, 2) +
                                   utils::powi(lscale_b, 2)))
                            : lscale_s *
                                  std::sqrt(
                                      1.0_rt -
                                      (utils::powi(Cmu, 6) *
                                       utils::powi(Cb_unstable, -2) * Rt));
                    tlscale_arrs[nbx](i, j, k) =
                        (stratification > 0)
                            ? amrex::min<amrex::Real>(
                                  tlscale_arrs[nbx](i, j, k),
                                  std::sqrt(
                                      Cmu * tke_arrs[nbx](i, j, k) /
                                      stratification))
                            : tlscale_arrs[nbx](i, j, k);
                    tlscale_arrs[nbx](i, j, k) =
                        (std::abs(surf_flux) < 1.0e-5_rt &&
                         z <= lengthscale_switch)
                            ? lscale_s
                            : tlscale_arrs[nbx](i, j, k);
                    Rt = (std::abs(surf_flux) < 1.0e-5_rt &&
                          z <= lengthscale_switch)
                             ? 0.0_rt
                             : Rt;
                    const amrex::Real Cmu_Rt =
                        (Cmu + (0.108_rt * Rt)) /
                        (1.0_rt + (0.308_rt * Rt) +
                         (0.00837_rt * utils::powi(Rt, 2)));
                    mu_arrs[nbx](i, j, k) = rho_arrs[nbx](i, j, k) * Cmu_Rt *
                                            tlscale_arrs[nbx](i, j, k) *
                                            std::sqrt(tke_arrs[nbx](i, j, k));
                    const amrex::Real Cmu_prime_Rt =
                        Cmu / (1.0_rt + (0.277_rt * Rt));
                    const amrex::Real muPrime =
                        rho_arrs[nbx](i, j, k) * Cmu_prime_Rt *
                        tlscale_arrs[nbx](i, j, k) *
                        std::sqrt(tke_arrs[nbx](i, j, k));
                    buoy_prod_arrs[nbx](i, j, k) = -muPrime * stratification;
                    shear_prod_arrs[nbx](i, j, k) *=
                        shear_prod_arrs[nbx](i, j, k) * mu_arrs[nbx](i, j, k);
                });
        }
    }
    amrex::Gpu::streamSynchronize();

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

// TerrainDrag, KLAxell.terrain_wall_stencil: the central strain rate of a
//  fluid cell next to a blanked cell reads the near-zero velocity the drag
//  holds in the blanked cell, not a fluid value. For a wall on the face
//  between them dU/dz = U_{k+1} / (2 dz), 1.4 times the log-law shear
//  u* / (kappa dz / 2), and the wall-cell TKE settles 16-22 % above the
//  target KransAxell relaxes it to, where flat ground sits 9-12 % below it.
//  Flat ground takes a one-sided stencil at its bottom wall (StencilKLO with
//  the hoextrap ghost of the tangential velocity and the zero wall value of
//  the normal one):
//      du_t/dn = (-3 u_0 + 4 u_1 - u_2) / (2 dx),
//      du_n/dn = (u_0 + u_1 / 3) / dx,
//  with u_0 the wall cell and u_1, u_2 the next fluid cells away from the
//  wall. The cells next to the blanked cells take the same stencil in each
//  direction that has a blanked neighbor, so that a wall on a cell face
//  reproduces flat ground, error included.
template <typename Transport>
void KLAxell<Transport>::terrain_wall_strain_rate(const Field& velocity)
{
    const auto& repo = this->m_sim.repo();
    const auto& blank = repo.get_int_field("terrain_blank");
    auto& strain = this->m_shear_prod;
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = repo.mesh().Geom(lev);
        const auto dx = geom.CellSizeArray();
        const auto dlo = amrex::lbound(geom.Domain());
        const auto dhi = amrex::ubound(geom.Domain());
        const amrex::GpuArray<int, AMREX_SPACEDIM> periodic{
            static_cast<int>(geom.isPeriodic(0)),
            static_cast<int>(geom.isPeriodic(1)),
            static_cast<int>(geom.isPeriodic(2))};
        const auto& s_arrs = strain(lev).arrays();
        const auto& v_arrs = velocity(lev).const_arrays();
        const auto& b_arrs = blank(lev).const_arrays();
        amrex::ParallelFor(
            strain(lev),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const auto& bl = b_arrs[nbx];
                if (bl(i, j, k) == 1) {
                    return;
                }
                const amrex::GpuArray<int, AMREX_SPACEDIM> idx{i, j, k};
                const amrex::GpuArray<int, AMREX_SPACEDIM> lo{
                    dlo.x, dlo.y, dlo.z};
                const amrex::GpuArray<int, AMREX_SPACEDIM> hi{
                    dhi.x, dhi.y, dhi.z};
                amrex::GpuArray<int, 2 * AMREX_SPACEDIM> side{};
                bool wall = false;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    const int di = (d == 0) ? 1 : 0;
                    const int dj = (d == 1) ? 1 : 0;
                    const int dk = (d == 2) ? 1 : 0;
                    side[2 * d] = ((periodic[d] == 0) && (idx[d] == lo[d]))
                                      ? 2
                                      : bl(i - di, j - dj, k - dk);
                    side[(2 * d) + 1] =
                        ((periodic[d] == 0) && (idx[d] == hi[d]))
                            ? 2
                            : bl(i + di, j + dj, k + dk);
                    wall =
                        wall || (side[2 * d] == 1) || (side[(2 * d) + 1] == 1);
                }
                if (!wall) {
                    return;
                }
                s_arrs[nbx](i, j, k) =
                    wall_strain_rate(v_arrs[nbx], i, j, k, side, dx);
            });
    }
    amrex::Gpu::streamSynchronize();
}

// TerrainDrag, KLAxell.terrain_blanked_face_length: a cell is blanked when
//  its center is at or below the terrain height h, so the resolved wall is
//  the top face of the blanked column,
//      z_f = z_lo + max(0, dz floor((h - z_lo) / dz + 1/2)),
//  not h. The mixing length of the TerrainDrag kernel is measured from z_f,
//  z = max(z_c - z_f, dz / 2), as flat ground measures it from its wall; with
//  h between cell faces the height above h puts the first fluid cells up to
//  dz / 2 too close to the wall.
template <typename Transport>
void KLAxell<Transport>::terrain_blanked_face_height(
    ScratchField& face_height) const
{
    const auto& repo = this->m_sim.repo();
    const auto& terrain_height = repo.get_field("terrain_height");
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = repo.mesh().Geom(lev);
        const amrex::Real zlo = geom.ProbLoArray()[2];
        const amrex::Real dz = geom.CellSize()[2];
        const auto& h_arrs = terrain_height(lev).const_arrays();
        const auto& f_arrs = face_height(lev).arrays();
        amrex::ParallelFor(
            face_height(lev),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
                const amrex::Real h = h_arrs[nbx](i, j, k) - zlo;
                f_arrs[nbx](i, j, k) =
                    zlo + amrex::max<amrex::Real>(
                              0.0_rt, dz * std::floor((h / dz) + 0.5_rt));
            });
    }
    amrex::Gpu::streamSynchronize();
}

// TerrainDrag, KLAxell.terrain_face_stress: DragForcing relaxes the drag
//  cell (terrain_drag == 1, first fluid cell above the blanked column) toward
//  the log-law velocity of the cell above it, which fixes the velocity
//  difference across the face between the two cells; the resolved flux
//  through that face is then the wall stress. The blanked cell below has no
//  turbulent viscosity, so the model viscosity leaves that flux short of
//  u*^2. Here the drag-cell viscosity is sized so that the face viscosity,
//  the mean of the two cells, carries the wall stress of the DragForcing
//  wall law over the current velocity difference:
//      u* = kappa |U_{k+1}| / (ln(1.5 dz / z0) - psi_m(1.5 dz / L)),
//      mu_f = rho u*^2 dz / max(|U_{k+1} - U_k|, 1e-2),
//      mu_k = max(2 mu_f - mu_{k+1}, 0),
//  with z0 floored at DragForcing.minimum_z0 and psi_m = 0 unless
//  ABL.wall_het_model = mol. The shear production of the drag cell keeps the
//  model viscosity.
template <typename Transport>
void KLAxell<Transport>::terrain_face_stress(
    const int lev, const Field& velocity, const Field& density)
{
    const auto& repo = this->m_sim.repo();
    const amrex::Real dz = repo.mesh().Geom(lev).CellSize()[2];
    const amrex::Real kappa = m_kappa;
    const amrex::Real min_z0 = m_min_z0;
    const amrex::Real psi_ref =
        (m_wall_het_model == "mol")
            ? MOData::calc_psi_m(
                  1.5_rt * dz / m_monin_obukhov_length, m_beta_m, m_gamma_m)
            : 0.0_rt;
    auto& mu_turb = this->mu_turb();
    const auto& mu_arrs = mu_turb(lev).arrays();
    const auto& rho_arrs = density(lev).const_arrays();
    const auto& v_arrs = velocity(lev).const_arrays();
    const auto& drag_arrs =
        repo.get_int_field("terrain_drag")(lev).const_arrays();
    const auto& z0_arrs = repo.get_field("terrainz0")(lev).const_arrays();
    amrex::ParallelFor(
        mu_turb(lev),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            if (drag_arrs[nbx](i, j, k) != 1) {
                return;
            }
            const auto& v = v_arrs[nbx];
            const amrex::Real z0 =
                amrex::max<amrex::Real>(z0_arrs[nbx](i, j, k), min_z0);
            const amrex::Real mref = std::sqrt(
                (v(i, j, k + 1, 0) * v(i, j, k + 1, 0)) +
                (v(i, j, k + 1, 1) * v(i, j, k + 1, 1)));
            const amrex::Real us =
                kappa * mref / (std::log(1.5_rt * dz / z0) - psi_ref);
            const amrex::Real du = v(i, j, k + 1, 0) - v(i, j, k, 0);
            const amrex::Real dv = v(i, j, k + 1, 1) - v(i, j, k, 1);
            const amrex::Real dM = amrex::max<amrex::Real>(
                std::sqrt((du * du) + (dv * dv)), 1.0e-2_rt);
            const amrex::Real rho = rho_arrs[nbx](i, j, k);
            const amrex::Real mu_face = rho * us * us * dz / dM;
            mu_arrs[nbx](i, j, k) = amrex::max<amrex::Real>(
                (2.0_rt * mu_face) - mu_arrs[nbx](i, j, k + 1), 0.0_rt);
        });
}

template <typename Transport>
void KLAxell<Transport>::update_alphaeff(Field& alphaeff)
{

    BL_PROFILE("kynema-sgf::" + this->identifier() + "::update_alphaeff");
    auto lam_alpha = (this->m_transport).alpha();
    auto& mu_turb = this->m_mu_turb;
    auto& repo = mu_turb.repo();

    fvm::gradient(*m_gradT, m_temperature);
    auto& gradT = *m_gradT;
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};
    const auto beta = (this->m_transport).beta();
    const amrex::Real Cmu = m_Cmu;
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& muturb_arrs = mu_turb(lev).arrays();
        const auto& alphaeff_arrs = alphaeff(lev).arrays();
        const auto& lam_diff_arrs = (*lam_alpha)(lev).arrays();
        const auto& tke_arrs = (*this->m_tke)(lev).arrays();
        const auto& gradT_arrs = gradT(lev).const_arrays();
        const auto& tlscale_arrs = (this->m_turb_lscale)(lev).arrays();
        const auto& beta_arrs = (*beta)(lev).const_arrays();
        const amrex::Real Rtc = -1.0_rt;
        const amrex::Real Rtmin = -3.0_rt;
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
                Rt = (Rt > Rtc) ? Rt
                                : amrex::max<amrex::Real>(
                                      Rt, Rt - (utils::powi(Rt - Rtc, 2) /
                                                (Rt + Rtmin - (2.0_rt * Rtc))));
                const amrex::Real prandtlRt =
                    (1.0_rt + (0.193_rt * Rt)) / (1.0_rt + (0.0302_rt * Rt));
                alphaeff_arrs[nbx](i, j, k) =
                    lam_diff_arrs[nbx](i, j, k) +
                    (muturb_arrs[nbx](i, j, k) / prandtlRt);
            });
        if (m_terrain_face_heat_flux && (m_wall_het_model == "mol") &&
            this->m_sim.repo().int_field_exists("terrain_drag")) {
            terrain_face_heat_flux(lev, alphaeff, *lam_alpha);
        }
    }
    amrex::Gpu::streamSynchronize();

    alphaeff.fillpatch(this->m_sim.time().current_time());
}

// TerrainDrag, KLAxell.terrain_face_heat_flux (ABL.wall_het_model = mol):
//  DragTempForcing relaxes the drag cell toward the MOST temperature from the
//  cell above it, which fixes the temperature difference across the face
//  between the two, so that face carries the surface heat flux. The drag-cell
//  heat diffusivity is sized so that the face (the mean of the two cells)
//  carries the MOST heat flux given by the Obukhov length, with the friction
//  velocity of terrain_face_stress:
//      theta* = theta_k u*^2 / (kappa g L),   q = -u* theta*,
//      alpha_f = rho q dz / (theta_k - theta_{k+1}),
//      alpha_k = max(2 alpha_f - alpha_{k+1}, alpha_lam).
//  Left unchanged where the difference is below 1e-4 K or runs against the
//  flux (q (theta_k - theta_{k+1}) <= 0).
template <typename Transport>
void KLAxell<Transport>::terrain_face_heat_flux(
    const int lev, Field& alphaeff, const ScratchField& lam_alpha)
{
    const auto& repo = this->m_sim.repo();
    const amrex::Real dz = repo.mesh().Geom(lev).CellSize()[2];
    const amrex::Real kappa = m_kappa;
    const amrex::Real min_z0 = m_min_z0;
    const amrex::Real L = m_monin_obukhov_length;
    const amrex::Real gmod = std::abs(m_gravity[2]);
    const amrex::Real psi_ref =
        MOData::calc_psi_m(1.5_rt * dz / L, m_beta_m, m_gamma_m);
    const auto& alpha_arrs = alphaeff(lev).arrays();
    const auto& lam_arrs = lam_alpha(lev).const_arrays();
    const auto& v_arrs = m_vel(lev).const_arrays();
    const auto& t_arrs = m_temperature(lev).const_arrays();
    const auto& rho_arrs = m_rho(lev).const_arrays();
    const auto& drag_arrs =
        repo.get_int_field("terrain_drag")(lev).const_arrays();
    const auto& z0_arrs = repo.get_field("terrainz0")(lev).const_arrays();
    amrex::ParallelFor(
        alphaeff(lev),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            if (drag_arrs[nbx](i, j, k) != 1) {
                return;
            }
            const auto& v = v_arrs[nbx];
            const auto& t = t_arrs[nbx];
            const amrex::Real z0 =
                amrex::max<amrex::Real>(z0_arrs[nbx](i, j, k), min_z0);
            const amrex::Real mref = std::sqrt(
                (v(i, j, k + 1, 0) * v(i, j, k + 1, 0)) +
                (v(i, j, k + 1, 1) * v(i, j, k + 1, 1)));
            const amrex::Real us =
                kappa * mref / (std::log(1.5_rt * dz / z0) - psi_ref);
            const amrex::Real ts = t(i, j, k) * us * us / (kappa * gmod * L);
            // Upward heat flux, positive when the surface heats the air
            const amrex::Real q = -us * ts;
            // Down-gradient: the flux runs from the drag cell up
            const amrex::Real dth = t(i, j, k) - t(i, j, k + 1);
            if ((q * dth <= 0.0_rt) || (std::abs(dth) < 1.0e-4_rt)) {
                return;
            }
            const amrex::Real a_face = rho_arrs[nbx](i, j, k) * q * dz / dth;
            alpha_arrs[nbx](i, j, k) = amrex::max<amrex::Real>(
                (2.0_rt * a_face) - alpha_arrs[nbx](i, j, k + 1),
                lam_arrs[nbx](i, j, k));
        });
}

template <typename Transport>
void KLAxell<Transport>::update_scalar_diff(
    Field& deff, const std::string& name)
{
    BL_PROFILE("kynema-sgf::" + this->identifier() + "::update_scalar_diff");

    if (name == pde::TKE::var_name()) {
        auto& mu_turb = this->mu_turb();
        deff.setVal(0.0_rt);
        field_ops::saxpy(
            deff, 2.0_rt, mu_turb, 0, 0, deff.num_comp(), deff.num_grow());
    } else {
        amrex::Abort(
            "KLAxell:update_scalar_diff not implemented for field " + name);
    }
}

template <typename Transport>
void KLAxell<Transport>::post_advance_work()
{
    BL_PROFILE("kynema-sgf::" + this->identifier() + "::post_advance_work");
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(KLAxell);

} // namespace kynema_sgf
