#include "amr-wind/turbulence/RANS/KOmegaSST.H"
#include "amr-wind/turbulence/RANS/KOmegaSSTI.H"
#include "amr-wind/equation_systems/PDEBase.H"
#include "amr-wind/turbulence/TurbModelDefs.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/fvm/strainrate.H"
#include "amr-wind/turbulence/turb_utils.H"
#include "amr-wind/equation_systems/tke/TKE.H"
#include "amr-wind/equation_systems/sdr/SDR.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {
namespace turbulence {

template <typename Transport>
void KOmegaSST<Transport>::parse_model_coeffs()
{

    if ((!this->m_sim.pde_manager().constant_density() ||
         this->m_sim.physics_manager().contains("MultiPhase"))) {
        this->m_include_buoyancy = true;
    }
    this->m_buoyancy_factor = (this->m_include_buoyancy) ? 1.0_rt : 0.0_rt;

    {
        const std::string coeffs_dict = this->model_name() + "_coeffs";
        amrex::ParmParse pp(coeffs_dict);
        pp.query("beta_star", this->m_beta_star);
        pp.query("alpha1", this->m_alpha1);
        pp.query("alpha2", this->m_alpha2);
        pp.query("beta1", this->m_beta1);
        pp.query("beta2", this->m_beta2);
        pp.query("sigma_k1", this->m_sigma_k1);
        pp.query("sigma_k2", this->m_sigma_k2);
        pp.query("sigma_omega1", this->m_sigma_omega1);
        pp.query("sigma_omega2", this->m_sigma_omega2);
        pp.query("tke_amb", this->m_tke_amb);
        pp.query("sdr_amb", this->m_sdr_amb);
        pp.query("sigma_t", this->m_sigma_t);
    }

    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }
}

template <typename Transport>
TurbulenceModel::CoeffsDictType KOmegaSST<Transport>::model_coeffs() const
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
        {"a1", this->m_a1}};
}

template <typename Transport>
void KOmegaSST<Transport>::update_turbulent_viscosity(
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

    auto& mu_turb = this->mu_turb();
    auto lam_mu = (this->m_transport).mu();
    const auto& den = this->m_rho.state(fstate);
    const auto& tke = (*this->m_tke).state(fstate);
    const auto& sdr = (*this->m_sdr).state(fstate);
    const auto& repo = mu_turb.repo();
    auto& tke_lhs = (this->m_sim).repo().get_field("tke_lhs_src_term");
    tke_lhs.setVal(0.0_rt);
    auto& sdr_lhs = (this->m_sim).repo().get_field("sdr_lhs_src_term");

    auto gradK = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradK, tke);

    auto gradOmega = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradOmega, sdr);

    // This is used for the buoyancy-modified version of the model
    auto gradden = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradden, den);

    const auto& vel = this->m_vel.state(fstate);
    // Compute strain rate into shear production term
    fvm::strainrate(this->m_shear_prod, vel);

    const amrex::Real delta_t = (this->m_sim).time().delta_t();
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};
    const amrex::Real Bfac = this->m_buoyancy_factor;
    const amrex::Real sigmat = this->m_sigma_t;

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& lam_mu_arrs = (*lam_mu)(lev).const_arrays();
        const auto& mu_arrs = mu_turb(lev).arrays();
        const auto& gradrho_arrs = (*gradden)(lev).const_arrays();
        const auto& rho_arrs = den(lev).const_arrays();
        const auto& gradK_arrs = (*gradK)(lev).const_arrays();
        const auto& gradOmega_arrs = (*gradOmega)(lev).const_arrays();
        const auto& tke_arrs = tke(lev).const_arrays();
        const auto& sdr_arrs = sdr(lev).const_arrays();
        const auto& wd_arrs = (this->m_walldist)(lev).const_arrays();
        const auto& shear_prod_arrs = (this->m_shear_prod)(lev).arrays();
        const auto& diss_arrs = (this->m_diss)(lev).arrays();
        const auto& sdr_src_arrs = (this->m_sdr_src)(lev).arrays();
        const auto& sdr_diss_arrs = (this->m_sdr_diss)(lev).arrays();
        const auto& f1_arrs = (this->m_f1)(lev).arrays();
        const auto& tke_lhs_arrs = tke_lhs(lev).arrays();
        const auto& sdr_lhs_arrs = sdr_lhs(lev).arrays();
        const auto& buoy_arrs = (this->m_buoy_term(lev)).arrays();

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

                // Buoyancy term
                amrex::Real tmpB =
                    -((gravity[0] * gradrho_arrs[nbx](i, j, k, 0)) +
                      (gravity[1] * gradrho_arrs[nbx](i, j, k, 1)) +
                      (gravity[2] * gradrho_arrs[nbx](i, j, k, 2)));

                buoy_arrs[nbx](i, j, k) =
                    Bfac * tmpB *
                    (mu_arrs[nbx](i, j, k) / rho_arrs[nbx](i, j, k)) / sigmat;

                f1_arrs[nbx](i, j, k) = tmp_f1;

                // For TKE equation:
                shear_prod_arrs[nbx](i, j, k) = amrex::min<amrex::Real>(
                    amrex::max<amrex::Real>(
                        mu_arrs[nbx](i, j, k) * tmp4 * tmp4, 0.0_rt),
                    10.0_rt * beta_star * rho_arrs[nbx](i, j, k) *
                        tke_arrs[nbx](i, j, k) * sdr_arrs[nbx](i, j, k));

                const amrex::Real diss_amb =
                    beta_star * rho_arrs[nbx](i, j, k) * sdr_amb * tke_amb;

                diss_arrs[nbx](i, j, k) =
                    (-beta_star * rho_arrs[nbx](i, j, k) *
                     tke_arrs[nbx](i, j, k) * sdr_arrs[nbx](i, j, k)) +
                    diss_amb;

                tke_lhs_arrs[nbx](i, j, k) = beta_star *
                                             rho_arrs[nbx](i, j, k) *
                                             sdr_arrs[nbx](i, j, k) * delta_t;

                // For SDR equation:
                amrex::Real production_omega =
                    rho_arrs[nbx](i, j, k) * alpha *
                    amrex::min<amrex::Real>(
                        tmp4 * tmp4, 10.0_rt * beta_star *
                                         sdr_arrs[nbx](i, j, k) *
                                         sdr_arrs[nbx](i, j, k));

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

template <typename Transport>
void KOmegaSST<Transport>::update_scalar_diff(
    Field& deff, const std::string& name)
{
    BL_PROFILE("amr-wind::" + this->identifier() + "::update_scalar_diff");

    auto lam_mu = (this->m_transport).mu();
    const auto& mu_turb = this->mu_turb();

    if (name == pde::TKE::var_name()) {
        const amrex::Real sigma_k1 = this->m_sigma_k1;
        const amrex::Real sigma_k2 = this->m_sigma_k2;
        const auto& repo = deff.repo();
        const int nlevels = repo.num_active_levels();
        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& lam_mu_arrs = (*lam_mu)(lev).const_arrays();
            const auto& mu_arrs = mu_turb(lev).const_arrays();
            const auto& f1_arrs = (this->m_f1)(lev).const_arrays();
            const auto& deff_arrs = deff(lev).arrays();
            amrex::ParallelFor(
                deff(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    deff_arrs[nbx](i, j, k) =
                        lam_mu_arrs[nbx](i, j, k) +
                        ((f1_arrs[nbx](i, j, k) * (sigma_k1 - sigma_k2) +
                          sigma_k2) *
                         mu_arrs[nbx](i, j, k));
                });
        }
        amrex::Gpu::streamSynchronize();
    } else if (name == pde::SDR::var_name()) {
        const amrex::Real sigma_omega1 = this->m_sigma_omega1;
        const amrex::Real sigma_omega2 = this->m_sigma_omega2;
        const auto& repo = deff.repo();
        const int nlevels = repo.num_active_levels();
        for (int lev = 0; lev < nlevels; ++lev) {
            const auto& lam_mu_arrs = (*lam_mu)(lev).const_arrays();
            const auto& mu_arrs = mu_turb(lev).const_arrays();
            const auto& f1_arrs = (this->m_f1)(lev).const_arrays();
            const auto& deff_arrs = deff(lev).arrays();
            amrex::ParallelFor(
                deff(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                    deff_arrs[nbx](i, j, k) =
                        lam_mu_arrs[nbx](i, j, k) +
                        ((f1_arrs[nbx](i, j, k) *
                              (sigma_omega1 - sigma_omega2) +
                          sigma_omega2) *
                         mu_arrs[nbx](i, j, k));
                });
        }
        amrex::Gpu::streamSynchronize();
    } else {
        amrex::Abort(
            "KOmegaSST:update_scalar_diff not implemented for field " + name);
    }
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(KOmegaSST);

} // namespace amr_wind
