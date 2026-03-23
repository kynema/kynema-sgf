#ifndef AMD_H
#define AMD_H

#include <AMReX_GpuContainers.H>
#include <string>
#include "amr-wind/turbulence/TurbModelBase.H"
#include "amr-wind/utilities/FieldPlaneAveragingFine.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/fvm/stencils.H"
#include "amr-wind/utilities/DirectionSelector.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::turbulence {
/** AMD LES Model
 *  \ingroup turb_model
 */
template <typename Transport>
class AMD : public TurbModelBase<Transport>
{
public:
    static std::string identifier() { return "AMD-" + Transport::identifier(); }

    explicit AMD(CFDSim& sim);

    //! Model name for debugging purposes
    [[nodiscard]] std::string model_name() const override { return "AMD"; }

    //! Update the turbulent viscosity field
    void update_turbulent_viscosity(
        FieldState fstate, DiffusionType /*unused*/) override;

    template <typename IndexSelector>
    void
    update_turbulent_viscosity(FieldState fstate, const IndexSelector& idxOp);

    //! Update the effective thermal diffusivity field
    void update_alphaeff(Field& alphaeff) override;

    //! Return model coefficients dictionary
    [[nodiscard]] TurbulenceModel::CoeffsDictType model_coeffs() const override;

    //! Parse turbulence model coefficients for this model
    void parse_model_coeffs() override;

    //! No post advance work for this model
    void post_advance_work() override {};

private:
    //! Poincare coefficient (default value set for 2nd order AMR-wind
    //! discretization)
    amrex::Real m_C{0.333333333333333_rt};

    //! Wall-normal direction axis
    int m_normal_dir{2};

    const Field& m_vel;
    const Field& m_temperature;
    const Field& m_rho;
    FieldPlaneAveragingFine m_pa_temp;
    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};
};

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real amd_muvel(
    int i,
    int j,
    int k,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx, // Grid spacing
    const amrex::Real beta, // Thermal expansion coefficient
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& gravity,
    const amrex::Real C, // Poincare const
    const amrex::Array4<amrex::Real const>& gradVel,
    const amrex::Array4<amrex::Real const>& gradT,
    const amrex::Real* gradTbar_coord_begin,
    const amrex::Real* gradTbar_coord_end,
    const amrex::Real* gradTbar,
    const int normal_dir,
    const amrex::Real nlo)
{
    amrex::Real num_shear = 0.0_rt;
    amrex::Real num_buoy = 0.0_rt;
    amrex::Real denom = 0.0_rt;
    const int n_ind = amr_wind::direction_selector(i, j, k, normal_dir);
    const amrex::Real h = nlo + ((n_ind + 0.5_rt) * dx[normal_dir]);
    const amrex::Real gradTbar_h = amr_wind::interp::linear(
        gradTbar_coord_begin, gradTbar_coord_end, gradTbar, h);
    for (int ii = 0; ii < AMREX_SPACEDIM; ++ii) {
        // This should operate only on the wall normal velocity
        num_buoy +=
            gradVel(i, j, k, (normal_dir * AMREX_SPACEDIM) + ii) *
            (gradT(i, j, k, ii) - ((ii == normal_dir) ? gradTbar_h : 0.0_rt)) *
            dx[ii] * dx[ii];

        for (int jj = 0; jj < AMREX_SPACEDIM; ++jj) {
            const amrex::Real diuj =
                gradVel(i, j, k, (ii * AMREX_SPACEDIM) + jj);
            const amrex::Real djui =
                gradVel(i, j, k, (jj * AMREX_SPACEDIM) + ii);
            denom += diuj * diuj;
            const amrex::Real sij = 0.5_rt * (diuj + djui);
            for (int kk = 0; kk < AMREX_SPACEDIM; ++kk) {
                const amrex::Real dkui =
                    gradVel(i, j, k, (ii * AMREX_SPACEDIM) + kk);
                const amrex::Real dkuj =
                    gradVel(i, j, k, (jj * AMREX_SPACEDIM) + kk);
                num_shear += dkui * dkuj * dx[kk] * dx[kk] * sij;
            }
        }
    }
    denom = amrex::max<amrex::Real>(
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e1_rt, denom);
    num_shear *= -C;
    num_buoy *= C * beta * std::abs(gravity[normal_dir]);
    return amrex::max<amrex::Real>(
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e1_rt,
        (num_shear + num_buoy) / denom);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real amd_thermal_diff(
    int i,
    int j,
    int k,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx, // Grid spacing
    const amrex::Real C,                                    // Poincare const
    const amrex::Array4<amrex::Real const>& gradVel,
    const amrex::Array4<amrex::Real const>& gradT)
{
    amrex::Real num = 0;
    amrex::Real denom = 0;
    for (int ii = 0; ii < AMREX_SPACEDIM; ++ii) {
        const amrex::Real diT = gradT(i, j, k, ii);
        denom += diT * diT;
        for (int kk = 0; kk < AMREX_SPACEDIM; ++kk) {
            const amrex::Real dkui =
                gradVel(i, j, k, (ii * AMREX_SPACEDIM) + kk);
            const amrex::Real dkT = gradT(i, j, k, kk);
            num += dkui * diT * dkT * dx[kk] * dx[kk];
        }
    }

    denom = amrex::max<amrex::Real>(
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e1_rt, denom);
    num *= -C;
    return amrex::max<amrex::Real>(
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e1_rt, num / denom);
}

} // namespace amr_wind::turbulence

#endif /* AMD_H */
