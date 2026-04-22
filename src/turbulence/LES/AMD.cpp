#include <AMReX_GpuContainers.H>
#include <AMReX_GpuQualifiers.H>
#include <cmath>

#include "amr-wind/fvm/gradient.H"
#include "amr-wind/turbulence/LES/AMD.H"
#include "amr-wind/turbulence/TurbModelDefs.H"
#include "amr-wind/utilities/DirectionSelector.H"

#include "AMReX_REAL.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParmParse.H"

using namespace amrex::literals;

namespace amr_wind {
namespace turbulence {

template <typename Transport>
// cppcheck-suppress uninitMemberVar
AMD<Transport>::AMD(CFDSim& sim)
    : TurbModelBase<Transport>(sim)
    , m_vel(sim.repo().get_field("velocity"))
    , m_temperature(sim.repo().get_field("temperature"))
    , m_rho(sim.repo().get_field("density"))
    , m_pa_temp(m_temperature, sim.time(), m_normal_dir, true)
{
    {
        amrex::ParmParse pp("incflo");
        pp.queryarr("gravity", m_gravity);
    }
}

template <typename Transport>
void AMD<Transport>::parse_model_coeffs()
{
    const std::string coeffs_dict = this->model_name() + "_coeffs";
    amrex::ParmParse pp(coeffs_dict);
    pp.query("C_poincare", m_C);
}

template <typename Transport>
void AMD<Transport>::update_turbulent_viscosity(
    const FieldState fstate, const DiffusionType /*unused*/)
{
    BL_PROFILE(
        "amr-wind::" + this->identifier() + "::update_turbulent_viscosity");

    auto& mu_turb = this->mu_turb();
    const auto& repo = mu_turb.repo();
    const auto& vel = m_vel.state(fstate);
    const auto& temp = m_temperature.state(fstate);
    const auto& den = m_rho.state(fstate);
    const auto beta = (this->m_transport).beta();
    const auto& geom_vec = repo.mesh().Geom();
    const int normal_dir = m_normal_dir;

    const amrex::Real C_poincare = m_C;

    auto gradVel = repo.create_scratch_field(AMREX_SPACEDIM * AMREX_SPACEDIM);
    fvm::gradient(*gradVel, vel);
    auto gradT = repo.create_scratch_field(AMREX_SPACEDIM);
    fvm::gradient(*gradT, temp);
    m_pa_temp(); // compute the current plane average
    const auto& tpa_deriv = m_pa_temp.line_deriv();
    amrex::Vector<amrex::Real> tpa_coord(tpa_deriv.size(), 0.0_rt);
    for (int i = 0; i < m_pa_temp.ncell_line(); ++i) {
        tpa_coord[i] = m_pa_temp.xlo() + ((0.5_rt + i) * m_pa_temp.dx());
    }
    amrex::Gpu::DeviceVector<amrex::Real> tpa_deriv_d(tpa_deriv.size(), 0.0_rt);
    amrex::Gpu::DeviceVector<amrex::Real> tpa_coord_d(tpa_coord.size(), 0.0_rt);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tpa_coord.begin(), tpa_coord.end(),
        tpa_coord_d.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, tpa_deriv.begin(), tpa_deriv.end(),
        tpa_deriv_d.begin());

    const amrex::Real* p_tpa_coord_begin = tpa_coord_d.data();
    const amrex::Real* p_tpa_coord_end = tpa_coord_d.end();
    const amrex::Real* p_tpa_deriv = tpa_deriv_d.data();
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> gravity{
        m_gravity[0], m_gravity[1], m_gravity[2]};
    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];
        const auto& problo = geom.ProbLoArray();
        const amrex::Real nlo = problo[normal_dir];
        const auto& dx = geom.CellSizeArray();

        const auto& gradVel_arrs = (*gradVel)(lev).const_arrays();
        const auto& gradT_arrs = (*gradT)(lev).const_arrays();
        const auto& rho_arrs = den(lev).const_arrays();
        const auto& beta_arrs = (*beta)(lev).const_arrays();
        const auto& mu_arrs = mu_turb(lev).arrays();
        const auto& mu_turb_lev = mu_turb(lev);
        amrex::ParallelFor(
            mu_turb_lev, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                auto mu_arr = mu_arrs[nbx];
                const auto rho_arr = rho_arrs[nbx];
                const auto gradVel_arr = gradVel_arrs[nbx];
                const auto gradT_arr = gradT_arrs[nbx];
                const auto beta_val = beta_arrs[nbx](i, j, k);
                mu_arr(i, j, k) =
                    rho_arr(i, j, k) * amd_muvel(
                                           i, j, k, dx, beta_val, gravity,
                                           C_poincare, gradVel_arr, gradT_arr,
                                           p_tpa_coord_begin, p_tpa_coord_end,
                                           p_tpa_deriv, normal_dir, nlo);
            });
    }
    amrex::Gpu::streamSynchronize();

    mu_turb.fillpatch(this->m_sim.time().current_time());
}

//! Update the effective thermal diffusivity field
template <typename Transport>
void AMD<Transport>::update_alphaeff(Field& alphaeff)
{

    BL_PROFILE("amr-wind::" + this->identifier() + "::update_alphaeff");

    const auto& repo = alphaeff.repo();
    const auto& geom_vec = repo.mesh().Geom();
    const amrex::Real C_poincare = m_C;
    auto gradVel = repo.create_scratch_field(AMREX_SPACEDIM * AMREX_SPACEDIM);
    fvm::gradient(*gradVel, m_vel);
    auto gradT = repo.create_scratch_field(AMREX_SPACEDIM);
    fvm::gradient(*gradT, m_temperature);

    const int nlevels = repo.num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        const auto& geom = geom_vec[lev];

        const auto& dx = geom.CellSizeArray();
        const auto& gradVel_arrs = (*gradVel)(lev).const_arrays();
        const auto& gradT_arrs = (*gradT)(lev).const_arrays();
        const auto& rho_arrs = m_rho(lev).const_arrays();
        const auto& alpha_arrs = alphaeff(lev).arrays();
        amrex::ParallelFor(
            alphaeff(lev), [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real rho = rho_arrs[nbx](i, j, k);
                alpha_arrs[nbx](i, j, k) =
                    rho * amd_thermal_diff(
                              i, j, k, dx, C_poincare, gradVel_arrs[nbx],
                              gradT_arrs[nbx]);
            });
    }
    amrex::Gpu::streamSynchronize();
}

template <typename Transport>
TurbulenceModel::CoeffsDictType AMD<Transport>::model_coeffs() const
{
    return TurbulenceModel::CoeffsDictType{{"C_poincare", m_C}};
}

} // namespace turbulence

INSTANTIATE_TURBULENCE_MODEL(AMD);

} // namespace amr_wind
