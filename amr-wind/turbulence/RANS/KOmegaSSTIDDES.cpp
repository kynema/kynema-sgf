#include "amr-wind/turbulence/RANS/KOmegaSSTIDDES.H"
#include "amr-wind/turbulence/RANS/KOmegaSSTI.H"
#include "amr-wind/equation_systems/PDEBase.H"
#include "amr-wind/turbulence/TurbModelDefs.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/fvm/strainrate.H"
#include "amr-wind/fvm/vorticity.H"
#include "amr-wind/fvm/vorticity_mag.H"
#include "amr-wind/turbulence/turb_utils.H"
#include "amr-wind/equation_systems/tke/TKE.H"
#include "amr-wind/equation_systems/sdr/SDR.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {
namespace turbulence {

template <typename Transport>
KOmegaSSTIDDES<Transport>::~KOmegaSSTIDDES() = default;

template <typename Transport>
// cppcheck-suppress uninitMemberVar
KOmegaSSTIDDES<Transport>::KOmegaSSTIDDES(CFDSim& sim)
    : KOmegaSST<Transport>(sim)
    , m_rans_ind(sim.repo().declare_field("rans_indicator", 1, 1, 1))
{}

template <typename Transport>
void KOmegaSSTIDDES<Transport>::parse_model_coeffs()
{
    KOmegaSST<Transport>::parse_model_coeffs();
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Cdes1", this->m_Cdes1);
    pp.query("Cdes2", this->m_Cdes2);
    pp.query("Cdt1", this->m_Cdt1);
    pp.query("Cdt2", this->m_Cdt2);
    pp.query("Cl", this->m_Cl);
    pp.query("Ct", this->m_Ct);
    pp.query("Cw", this->m_Cw);
    pp.query("kappa", this->m_kappa);
    pp.query("sdr_prod_clip_factor", this->m_sdr_prod_clip_factor);
}

template <typename Transport>
TurbulenceModel::CoeffsDictType KOmegaSSTIDDES<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{
        {"beta_star", this->m_beta_star},
        {"alpha1", this->m_alpha1},
        {"alpha2", this->m_alpha2},
        {"beta1", this->m_beta1},
        {"beta2", this->m_beta2},
        {"sigma_k1", this->m_sigma_k1},
        {"sigma_k2", this->m_sigma_k2},
        {"sigma_omega1", this->m_sigma_omega1},
        {"sigma_omega2", this->m_sigma_omega2},
        {"tke_amb", this->m_tke_amb},
        {"sdr_amb", this->m_sdr_amb},
        {"a1", this->m_a1},
        {"Cdes1", this->m_Cdes1},
        {"Cdes2", this->m_Cdes2},
        {"Cdt1", this->m_Cdt1},
        {"Cdt2", this->m_Cdt2},
        {"Cl", this->m_Cl},
        {"Ct", this->m_Ct},
        {"Cw", this->m_Cw},
        {"kappa", this->m_kappa},
        {"sdr_prod_clip_factor", this->m_sdr_prod_clip_factor}};
}

template <typename Transport>
void KOmegaSSTIDDES<Transport>::update_turbulent_viscosity(
    const FieldState fstate, const DiffusionType diff_type)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");

    const amrex::Real beta_star = this->m_beta_star;
    const amrex::Real alpha1 = this->m_alpha1;
    const amrex::Real alpha2 = this->m_alpha2;
    const amrex::Real beta1 = this->m_beta1;
    const amrex::Real beta2 = this->m_beta2;
    const amrex::Real sigma_omega2 = this->m_sigma_omega2;
    const amrex::Real a1 = this->m_a1;
    const amrex::Real tke_amb = this->m_tke_amb;
    const amrex::Real sdr_amb = this->m_sdr_amb;
    const amrex::Real Cdes1 = this->m_Cdes1;
    const amrex::Real Cdes2 = this->m_Cdes2;
    const amrex::Real Cdt1 = this->m_Cdt1;
    const amrex::Real Cdt2 = this->m_Cdt2;
    const amrex::Real Cl = this->m_Cl;
    const amrex::Real Ct = this->m_Ct;
    const amrex::Real Cw = this->m_Cw;
    const amrex::Real kappa = this->m_kappa;
    const amrex::Real sdr_prod_clip_factor = this->m_sdr_prod_clip_factor;

    auto& mu_turb = this->mu_turb();
    auto lam_mu = (this->m_transport).mu();
    const auto& den = this->m_rho.state(fstate);
    const auto& tke = (*this->m_tke).state(fstate);
    const auto& sdr = (*this->m_sdr).state(fstate);
    const auto& repo = mu_turb.repo();
    const auto& geom_vec = repo.mesh().Geom();
    auto& tke_lhs = (this->m_sim).repo().get_field("tke_lhs_src_term");
    tke_lhs.setVal(0.0_rt);
    auto& sdr_lhs = (this->m_sim).repo().get_field("sdr_lhs_src_term");

    auto gradK = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradK, tke);

    auto gradOmega = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradOmega, sdr);

    const auto& vel = this->m_vel.state(fstate);
    // Compute strain rate into shear production term
    fvm::strainrate(this->m_shear_prod, vel);

    auto vortmag = (this->m_sim.repo()).create_scratch_field(1, 0);
    fvm::vorticity_mag(*vortmag, vel);

    const amrex::Real delta_t = (this->m_sim).time().delta_t();

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];
        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real hmax =
            amrex::max<amrex::Real>(amrex::max<amrex::Real>(dx, dy), dz);
        const auto& lam_mu_arrs = (*lam_mu)(lev).const_arrays();
        const auto& mu_arrs = mu_turb(lev).arrays();
        const auto& rho_arrs = den(lev).const_arrays();
        const auto& gradK_arrs = (*gradK)(lev).const_arrays();
        const auto& gradOmega_arrs = (*gradOmega)(lev).const_arrays();
        const auto& tke_arrs = tke(lev).const_arrays();
        const auto& sdr_arrs = sdr(lev).const_arrays();
        const auto& wd_arrs = (this->m_walldist)(lev).const_arrays();
        const auto& shear_prod_arrs = (this->m_shear_prod)(lev).arrays();
        const auto& vortmag_arrs = (*vortmag)(lev).const_arrays();
        const auto& rans_ind_arrs = (this->m_rans_ind)(lev).arrays();
        const auto& diss_arrs = (this->m_diss)(lev).arrays();
        const auto& sdr_src_arrs = (this->m_sdr_src)(lev).arrays();
        const auto& sdr_diss_arrs = (this->m_sdr_diss)(lev).arrays();
        const auto& f1_arrs = (this->m_f1)(lev).arrays();
        const auto& tke_lhs_arrs = tke_lhs(lev).arrays();
        const auto& sdr_lhs_arrs = sdr_lhs(lev).arrays();

        amrex::ParallelFor(
            mu_turb(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                amrex::Real gko =
                    ((gradK_arrs[nbx](i, j, k, 0) *
                      gradOmega_arrs[nbx](i, j, k, 0)) +
                     (gradK_arrs[nbx](i, j, k, 1) *
                      gradOmega_arrs[nbx](i, j, k, 1)) +
                     (gradK_arrs[nbx](i, j, k, 2) *
                      gradOmega_arrs[nbx](i, j, k, 2)));

                amrex::Real cdkomega = amrex::max<amrex::Real>(
                    std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt,
                    2.0_rt * rho_arrs[nbx](i, j, k) * sigma_omega2 * gko /
                        (sdr_arrs[nbx](i, j, k) +
                         std::numeric_limits<amrex::Real>::epsilon() *
                             1.0e1_rt));

                amrex::Real tmp1 =
                    4.0_rt * rho_arrs[nbx](i, j, k) * sigma_omega2 *
                    tke_arrs[nbx](i, j, k) /
                    (cdkomega * wd_arrs[nbx](i, j, k) * wd_arrs[nbx](i, j, k));
                amrex::Real tmp2 =
                    std::sqrt(tke_arrs[nbx](i, j, k)) /
                    (beta_star * sdr_arrs[nbx](i, j, k) *
                         wd_arrs[nbx](i, j, k) +
                     std::numeric_limits<amrex::Real>::epsilon() * 1.0e1_rt);
                amrex::Real tmp3 =
                    500.0_rt * lam_mu_arrs[nbx](i, j, k) /
                    (wd_arrs[nbx](i, j, k) * wd_arrs[nbx](i, j, k) *
                         sdr_arrs[nbx](i, j, k) * rho_arrs[nbx](i, j, k) +
                     std::numeric_limits<amrex::Real>::epsilon() * 1.0e1_rt);
                amrex::Real tmp4 = shear_prod_arrs[nbx](i, j, k);
                amrex::Real tmp5 = vortmag_arrs[nbx](i, j, k);

                amrex::Real arg1 = amrex::min<amrex::Real>(
                    amrex::max<amrex::Real>(tmp2, tmp3), tmp1);
                amrex::Real tmp_f1 = std::tanh(arg1 * arg1 * arg1 * arg1);

                amrex::Real alpha = (tmp_f1 * (alpha1 - alpha2)) + alpha2;
                amrex::Real beta = (tmp_f1 * (beta1 - beta2)) + beta2;

                amrex::Real arg2 = amrex::max<amrex::Real>(2.0_rt * tmp2, tmp3);
                amrex::Real f2 = std::tanh(arg2 * arg2);

                mu_arrs[nbx](i, j, k) =
                    rho_arrs[nbx](i, j, k) * a1 * tke_arrs[nbx](i, j, k) /
                    amrex::max<amrex::Real>(
                        a1 * sdr_arrs[nbx](i, j, k), tmp4 * f2);

                f1_arrs[nbx](i, j, k) = tmp_f1;

                const amrex::Real alpha_des =
                    0.25_rt - (wd_arrs[nbx](i, j, k) / hmax);
                const amrex::Real denom =
                    (rho_arrs[nbx](i, j, k) * kappa * kappa *
                     wd_arrs[nbx](i, j, k) * wd_arrs[nbx](i, j, k) *
                     std::sqrt(0.5_rt * (tmp4 * tmp4 + tmp5 * tmp5)));
                const amrex::Real rdl = lam_mu_arrs[nbx](i, j, k) / denom;
                const amrex::Real rdt = mu_arrs[nbx](i, j, k) / denom;
                const amrex::Real fl =
                    std::tanh(std::pow(Cl * Cl * rdl, 10.0_rt));
                const amrex::Real ft =
                    std::tanh(std::pow(Ct * Ct * rdt, 3.0_rt));
                const amrex::Real fe1 =
                    (alpha < 0) ? 2.0_rt * std::exp(-9.0_rt * alpha * alpha)
                                : 2.0_rt * std::exp(-11.09_rt * alpha * alpha);
                const amrex::Real fe2 =
                    1.0_rt - amrex::max<amrex::Real>(ft, fl);
                const amrex::Real fe =
                    fe2 * amrex::max<amrex::Real>((fe1 - 1.0_rt), 0.0_rt);
                const amrex::Real fb = amrex::min<amrex::Real>(
                    2.0_rt * std::exp(-9.0_rt * alpha_des * alpha_des), 1.0_rt);
                const amrex::Real fdt =
                    1.0_rt - std::tanh(std::pow(Cdt1 * rdt, Cdt2));
                const amrex::Real fdtilde =
                    amrex::max<amrex::Real>((1.0_rt - fdt), fb);
                const amrex::Real cdes = (tmp_f1 * (Cdes1 - Cdes2)) + Cdes2;
                const amrex::Real l_les =
                    cdes * amrex::min<amrex::Real>(
                               Cw * amrex::max<amrex::Real>(
                                        wd_arrs[nbx](i, j, k), hmax),
                               hmax);
                const amrex::Real l_rans = std::sqrt(tke_arrs[nbx](i, j, k)) /
                                           (beta_star * sdr_arrs[nbx](i, j, k));
                const amrex::Real rans_ind = fdtilde * (1.0_rt + fe);
                rans_ind_arrs[nbx](i, j, k) = rans_ind;
                const amrex::Real l_iddes = amrex::max<amrex::Real>(
                    std::numeric_limits<amrex::Real>::epsilon(),
                    (rans_ind * l_rans) + ((1.0_rt - fdtilde) * l_les));

                const amrex::Real sqrt_tke_amb = std::sqrt(tke_amb);
                const amrex::Real l_sst_amb =
                    sqrt_tke_amb /
                    (beta_star *
                     amrex::max<amrex::Real>(
                         std::numeric_limits<amrex::Real>::epsilon(), sdr_amb));
                const amrex::Real l_iddes_amb = amrex::max<amrex::Real>(
                    std::numeric_limits<amrex::Real>::epsilon(),
                    (rans_ind * l_sst_amb) + ((1.0_rt - fdtilde) * l_les));
                const amrex::Real diss_amb = rho_arrs[nbx](i, j, k) * tke_amb *
                                             sqrt_tke_amb / l_iddes_amb;

                // For TKE equation:
                shear_prod_arrs[nbx](i, j, k) = amrex::min<amrex::Real>(
                    amrex::max<amrex::Real>(
                        mu_arrs[nbx](i, j, k) * tmp4 * tmp4, 0.0_rt),
                    10.0_rt * rho_arrs[nbx](i, j, k) * tke_arrs[nbx](i, j, k) *
                        std::sqrt(tke_arrs[nbx](i, j, k)) / l_iddes);

                diss_arrs[nbx](i, j, k) = (-rho_arrs[nbx](i, j, k) *
                                           std::sqrt(tke_arrs[nbx](i, j, k)) *
                                           tke_arrs[nbx](i, j, k) / l_iddes) +
                                          diss_amb;

                tke_lhs_arrs[nbx](i, j, k) = rho_arrs[nbx](i, j, k) *
                                             std::sqrt(tke_arrs[nbx](i, j, k)) /
                                             l_iddes * delta_t;

                // For SDR equation:
                amrex::Real production_omega =
                    rho_arrs[nbx](i, j, k) * alpha *
                    amrex::min<amrex::Real>(
                        tmp4 * tmp4, sdr_prod_clip_factor *
                                         amrex::max<amrex::Real>(
                                             sdr_arrs[nbx](i, j, k), 0.0_rt) *
                                         std::sqrt(tke_arrs[nbx](i, j, k)) /
                                         l_iddes);

                amrex::Real cross_diffusion =
                    (1.0_rt - tmp_f1) * 2.0_rt * rho_arrs[nbx](i, j, k) *
                    sigma_omega2 * gko /
                    (sdr_arrs[nbx](i, j, k) +
                     std::numeric_limits<amrex::Real>::epsilon() * 1.0e1_rt);

                const amrex::Real sdr_diss_amb =
                    beta * rho_arrs[nbx](i, j, k) * sdr_amb * sdr_amb;

                if (diff_type == DiffusionType::Crank_Nicolson) {

                    tke_lhs_arrs[nbx](i, j, k) =
                        0.5_rt * tke_lhs_arrs[nbx](i, j, k);

                    sdr_src_arrs[nbx](i, j, k) = production_omega;

                    sdr_diss_arrs[nbx](i, j, k) = cross_diffusion;

                    sdr_lhs_arrs[nbx](i, j, k) =
                        (rho_arrs[nbx](i, j, k) * beta *
                             sdr_arrs[nbx](i, j, k) +
                         0.5_rt * std::abs(cross_diffusion) /
                             (sdr_arrs[nbx](i, j, k) +
                              std::numeric_limits<amrex::Real>::epsilon() *
                                  1.0e1_rt)) *
                        delta_t;

                } else if (diff_type == DiffusionType::Implicit) {
                    /* Source term linearization is based on Florian
                       Menter's (1993) AIAA paper */
                    diss_arrs[nbx](i, j, k) = 0.0_rt;

                    sdr_src_arrs[nbx](i, j, k) = production_omega;

                    sdr_diss_arrs[nbx](i, j, k) = 0.0_rt;

                    sdr_lhs_arrs[nbx](i, j, k) =
                        (2.0_rt * rho_arrs[nbx](i, j, k) * beta *
                             sdr_arrs[nbx](i, j, k) +
                         std::abs(cross_diffusion) /
                             (sdr_arrs[nbx](i, j, k) +
                              std::numeric_limits<amrex::Real>::epsilon() *
                                  1.0e1_rt)) *
                        delta_t;
                } else {
                    sdr_src_arrs[nbx](i, j, k) =
                        production_omega + cross_diffusion;

                    sdr_diss_arrs[nbx](i, j, k) =
                        (-rho_arrs[nbx](i, j, k) * beta *
                         sdr_arrs[nbx](i, j, k) * sdr_arrs[nbx](i, j, k)) +
                        sdr_diss_amb;

                    sdr_lhs_arrs[nbx](i, j, k) =
                        0.5_rt * rho_arrs[nbx](i, j, k) * beta *
                        sdr_arrs[nbx](i, j, k) * delta_t;
                }
            });
    }
    amrex::Gpu::streamSynchronize();

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(KOmegaSSTIDDES);

} // namespace amr_wind
