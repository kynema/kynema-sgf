#include "amr-wind/turbulence/LES/OneEqKsgs.H"
#include "amr-wind/equation_systems/PDEBase.H"
#include "amr-wind/turbulence/TurbModelDefs.H"
#include "amr-wind/fvm/gradient.H"
#include "amr-wind/fvm/strainrate.H"
#include "amr-wind/turbulence/turb_utils.H"
#include "amr-wind/equation_systems/tke/TKE.H"

#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {
namespace turbulence {

template <typename Transport>
// cppcheck-suppress uninitMemberVar
OneEqKsgs<Transport>::OneEqKsgs(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_turb_lscale(sim.repo().declare_field("turb_lscale", 1))
    , m_shear_prod(sim.repo().declare_field("shear_prod", 1))
    , m_buoy_prod(sim.repo().declare_field("buoy_prod", 1))
    , m_dissip(sim.repo().declare_field("dissipation", 1))
    , m_rho(sim.repo().get_field("density"))
{
    auto& tke_eqn =
        sim.pde_manager().register_transport_pde(pde::TKE::pde_name());
    m_tke = &(tke_eqn.fields().field);
}

template <typename Transport>
OneEqKsgs<Transport>::~OneEqKsgs() = default;

template <typename Transport>
OneEqKsgsM84<Transport>::OneEqKsgsM84(CFDSim& sim)
    : OneEqKsgs<Transport>(sim)
    , m_temperature(sim.repo().get_field("temperature"))
{

    auto& phy_mgr = this->m_sim.physics_manager();
    if (!phy_mgr.contains("ABL")) {
        amrex::Abort("OneEqKsgsM84 model only works with ABL physics");
    }

    {
        amrex::ParmParse pp("ABL");
        pp.query("enable_hybrid_rl_mode", m_hybrid_rl);
    }

    if (m_hybrid_rl) {
        m_sdr = &(sim.repo().declare_field(
            "sdr", 1, (*this->m_tke).num_grow()[0], 1));

        m_sdr->set_default_fillpatch_bc(sim.time());
    }

    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }

    // TKE source term to be added to PDE
    turb_utils::inject_turbulence_src_terms(
        pde::TKE::pde_name(), {"KsgsM84Src"});
}

template <typename Transport>
OneEqKsgsM84<Transport>::~OneEqKsgsM84() = default;

template <typename Transport>
void OneEqKsgsM84<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Ceps", this->m_Ceps);
    pp.query("Ce", this->m_Ce);
}

template <typename Transport>
TurbulenceModel::CoeffsDictType OneEqKsgsM84<Transport>::model_coeffs() const
{
    // clang-format off
    return TurbulenceModel::CoeffsDictType{
        {"Ce", this->m_Ce}, {"Ceps", this->m_Ceps}};
    // clang-format on
}

template <typename Transport>
void OneEqKsgsM84<Transport>::update_turbulent_viscosity(
    const FieldState fstate, const DiffusionType /*unused*/)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");

    auto gradT = (this->m_sim.repo()).create_scratch_field(3, 0);
    fvm::gradient(*gradT, m_temperature.state(fstate));

    const auto& vel = this->m_vel.state(fstate);
    // Compute strain rate into shear production term
    fvm::strainrate(this->m_shear_prod, vel);

    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};
    const auto beta = (this->m_transport).beta();

    auto& mu_turb = this->mu_turb();
    const amrex::Real Ce = this->m_Ce;
    const auto& den = this->m_rho.state(fstate);
    const auto& repo = mu_turb.repo();
    const auto& geom_vec = repo.mesh().Geom();

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real ds = std::cbrt(dx * dy * dz);
        const auto& mu_arrs = mu_turb(lev).arrays();
        const auto& rho_arrs = den(lev).const_arrays();
        const auto& gradT_arrs = (*gradT)(lev).const_arrays();
        const auto& tlscale_arrs = (this->m_turb_lscale)(lev).arrays();
        const auto& tke_arrs = (*this->m_tke)(lev).const_arrays();
        const auto& buoy_prod_arrs = (this->m_buoy_prod)(lev).arrays();
        const auto& shear_prod_arrs = (this->m_shear_prod)(lev).arrays();
        const auto& beta_arrs = (*beta)(lev).const_arrays();

        amrex::ParallelFor(
            mu_turb(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                amrex::Real stratification =
                    -((gradT_arrs[nbx](i, j, k, 0) * gravity[0]) +
                      (gradT_arrs[nbx](i, j, k, 1) * gravity[1]) +
                      (gradT_arrs[nbx](i, j, k, 2) * gravity[2])) *
                    beta_arrs[nbx](i, j, k);
                if (stratification >
                    std::numeric_limits<amrex::Real>::epsilon() * 1.0e6_rt) {
                    tlscale_arrs[nbx](i, j, k) = amrex::min<amrex::Real>(
                        ds,
                        0.76_rt *
                            std::sqrt(tke_arrs[nbx](i, j, k) / stratification));
                } else {
                    tlscale_arrs[nbx](i, j, k) = ds;
                }

                mu_arrs[nbx](i, j, k) = rho_arrs[nbx](i, j, k) * Ce *
                                        tlscale_arrs[nbx](i, j, k) *
                                        std::sqrt(tke_arrs[nbx](i, j, k));

                buoy_prod_arrs[nbx](i, j, k) =
                    -mu_arrs[nbx](i, j, k) *
                    (1.0_rt + 2.0_rt * tlscale_arrs[nbx](i, j, k) / ds) *
                    stratification;

                shear_prod_arrs[nbx](i, j, k) *=
                    shear_prod_arrs[nbx](i, j, k) * mu_arrs[nbx](i, j, k);
            });
    }
    amrex::Gpu::streamSynchronize();

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

template <typename Transport>
void OneEqKsgsM84<Transport>::update_alphaeff(Field& alphaeff)
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::update_alphaeff");

    auto lam_alpha = (this->m_transport).alpha();
    auto& mu_turb = this->m_mu_turb;
    auto& repo = mu_turb.repo();
    const auto& geom_vec = repo.mesh().Geom();

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real ds = std::cbrt(dx * dy * dz);

        const auto& muturb_arrs = mu_turb(lev).const_arrays();
        const auto& alphaeff_arrs = alphaeff(lev).arrays();
        const auto& tlscale_arrs = this->m_turb_lscale(lev).const_arrays();
        const auto& lam_diff_arrs = (*lam_alpha)(lev).const_arrays();

        amrex::ParallelFor(
            mu_turb(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                alphaeff_arrs[nbx](i, j, k) =
                    lam_diff_arrs[nbx](i, j, k) +
                    (muturb_arrs[nbx](i, j, k) *
                     (1.0_rt + 2.0_rt * tlscale_arrs[nbx](i, j, k) / ds));
            });
    }
    amrex::Gpu::streamSynchronize();

    alphaeff.fillpatch(this->m_sim.time().current_time());
}

template <typename Transport>
void OneEqKsgsM84<Transport>::update_scalar_diff(
    Field& deff, const std::string& name)
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::update_scalar_diff");

    if (name == pde::TKE::var_name()) {
        auto& mu_turb = this->mu_turb();
        deff.setVal(0.0_rt);
        field_ops::saxpy(
            deff, 2.0_rt, mu_turb, 0, 0, deff.num_comp(), deff.num_grow());
    } else {
        amrex::Abort(
            "OneEqKsgsM84:update_scalar_diff not implemented for field " +
            name);
    }
}

template <typename Transport>
void OneEqKsgsM84<Transport>::post_advance_work()
{

    if (!m_hybrid_rl) {
        return;
    }

    BL_PROFILE("amr-wind::" + this->identifier() + "::post_advance_work");

    // Update sdr field based on sfs ke

    auto& tke = *(this->m_tke);
    auto& sdr = *(this->m_sdr);
    const amrex::Real Ce = this->m_Ce;

    auto& repo = tke.repo();
    const auto& geom_vec = repo.mesh().Geom();
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real ds = std::cbrt(dx * dy * dz);

        const auto& tke_arrs = tke(lev).const_arrays();
        const auto& sdr_arrs = sdr(lev).arrays();

        amrex::ParallelFor(
            tke(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                sdr_arrs[nbx](i, j, k) =
                    std::sqrt(tke_arrs[nbx](i, j, k)) / (Ce * ds);
            });
    }
    amrex::Gpu::streamSynchronize();
}

template <typename Transport>
OneEqKsgsS94<Transport>::OneEqKsgsS94(CFDSim& sim) : OneEqKsgs<Transport>(sim)
{
    // TKE source term to be added to PDE
    turb_utils::inject_turbulence_src_terms(
        pde::TKE::pde_name(), {"KsgsS94Src"});
}

template <typename Transport>
OneEqKsgsS94<Transport>::~OneEqKsgsS94() = default;

template <typename Transport>
void OneEqKsgsS94<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Ceps", this->m_Ceps);
}

template <typename Transport>
TurbulenceModel::CoeffsDictType OneEqKsgsS94<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{{"Ceps", this->m_Ceps}};
}

template <typename Transport>
void OneEqKsgsS94<Transport>::update_turbulent_viscosity(
    const FieldState // fstate
    /*unused*/,
    const DiffusionType /*unused*/)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");
}

template <typename Transport>
void OneEqKsgsS94<Transport>::update_scalar_diff(
    Field& deff, const std::string& name)
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::update_scalar_diff");

    if (name == pde::TKE::var_name()) {
        auto& mu_turb = this->mu_turb();
        deff.setVal(0.0_rt);
        field_ops::saxpy(
            deff, 2.0_rt, mu_turb, 0, 0, deff.num_comp(), deff.num_grow());
    } else {
        amrex::Abort(
            "OneEqKsgsM84:update_scalar_diff not implemented for field " +
            name);
    }
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(OneEqKsgsM84);
INSTANTIATE_TURBULENCE_MODEL(OneEqKsgsS94);

} // namespace amr_wind
