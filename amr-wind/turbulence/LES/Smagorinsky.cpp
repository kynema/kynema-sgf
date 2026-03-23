#include <cmath>

#include "amr-wind/turbulence/LES/Smagorinsky.H"
#include "amr-wind/turbulence/TurbModelDefs.H"
#include "amr-wind/fvm/strainrate.H"
#include "AMReX_REAL.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace turbulence {

template <typename Transport>
// cppcheck-suppress uninitMemberVar
Smagorinsky<Transport>::Smagorinsky(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_rho(sim.repo().get_field("density"))
{}

template <typename Transport>
void Smagorinsky<Transport>::update_turbulent_viscosity(
    const FieldState fstate, const DiffusionType /*unused*/)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");

    auto& mu_turb = this->mu_turb();
    const auto& repo = mu_turb.repo();
    const auto& vel = m_vel.state(fstate);
    const auto& den = m_rho.state(fstate);
    const auto& geom_vec = repo.mesh().Geom();
    const amrex::Real Cs_sqr = this->m_Cs * this->m_Cs;

    // Populate strainrate into the turbulent viscosity arrays to avoid creating
    // a temporary buffer
    fvm::strainrate(mu_turb, vel);

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const amrex::Real dx = geom.CellSize()[0];
        const amrex::Real dy = geom.CellSize()[1];
        const amrex::Real dz = geom.CellSize()[2];
        const amrex::Real ds = std::cbrt(dx * dy * dz);
        const amrex::Real ds_sqr = ds * ds;
        const amrex::Real smag_factor = Cs_sqr * ds_sqr;

        const auto& mu_arrs = mu_turb(lev).arrays();
        const auto& rho_arrs = den(lev).const_arrays();

        amrex::ParallelFor(
            mu_turb(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real rho = rho_arrs[nbx](i, j, k);
                mu_arrs[nbx](i, j, k) *= rho * smag_factor;
            });
    }
    amrex::Gpu::streamSynchronize();

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

template <typename Transport>
void Smagorinsky<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("Cs", this->m_Cs);
}

template <typename Transport>
TurbulenceModel::CoeffsDictType Smagorinsky<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{{"Cs", this->m_Cs}};
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(Smagorinsky);

} // namespace amr_wind
