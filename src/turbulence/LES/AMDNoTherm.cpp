#include <cmath>

#include "amr-wind/fvm/gradient.H"
#include "amr-wind/turbulence/LES/AMDNoTherm.H"
#include "amr-wind/turbulence/TurbModelDefs.H"

#include "AMReX_Config.H"
#include "AMReX_REAL.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParmParse.H"

namespace amr_wind {
namespace turbulence {

template <typename Transport>
// cppcheck-suppress uninitMemberVar
AMDNoTherm<Transport>::AMDNoTherm(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_rho(sim.repo().get_field("density"))
{}

template <typename Transport>
void AMDNoTherm<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("C_poincare", m_C);
}

template <typename Transport>
void AMDNoTherm<Transport>::update_turbulent_viscosity(
    const FieldState fstate, const DiffusionType /*unused*/)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");

    auto& mu_turb = this->mu_turb();
    const auto& repo = mu_turb.repo();
    const auto& vel = m_vel.state(fstate);
    const auto& den = m_rho.state(fstate);
    auto gradVel = repo.create_scratch_field(AMREX_SPACEDIM * AMREX_SPACEDIM);
    fvm::gradient(*gradVel, vel);
    const auto& geom_vec = repo.mesh().Geom();

    const amrex::Real C_poincare = this->m_C;

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const auto& dx = geom.CellSizeArray();
        const auto& gradVel_arrs = (*gradVel)(lev).const_arrays();
        const auto& mu_arrs = mu_turb(lev).arrays();
        const auto& rho_arrs = den(lev).const_arrays();
        amrex::ParallelFor(
            mu_turb(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real rho = rho_arrs[nbx](i, j, k);
                mu_arrs[nbx](i, j, k) =
                    rho *
                    amd_base_muvel(i, j, k, dx, C_poincare, gradVel_arrs[nbx]);
            });
    }
    amrex::Gpu::streamSynchronize();

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

template <typename Transport>
TurbulenceModel::CoeffsDictType AMDNoTherm<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{{"C_poincare", this->m_C}};
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(AMDNoTherm);

} // namespace amr_wind
