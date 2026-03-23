#ifndef KOSOVIC_H
#define KOSOVIC_H

#include <string>
#include "amr-wind/turbulence/TurbModelBase.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::turbulence {

/** Kosovic LES Model
 *  \ingroup turb_model
 */
template <typename Transport>
class Kosovic : public TurbModelBase<Transport>
{
public:
    static std::string identifier()
    {
        return "Kosovic-" + Transport::identifier();
    }

    explicit Kosovic(CFDSim& sim);

    //! Model name for debugging purposes
    [[nodiscard]] std::string model_name() const override { return "Kosovic"; }

    //! Update the turbulent viscosity field
    void update_turbulent_viscosity(
        FieldState fstate, DiffusionType /*unused*/) override;

    //! No post advance work for this model
    void post_advance_work() override {}

    //! Update the effective thermal diffusivity field
    void update_alphaeff(Field& alphaeff) override;

    //! Parse turbulence model coefficients
    void parse_model_coeffs() override;

    //! Return model coefficients dictionary
    [[nodiscard]] TurbulenceModel::CoeffsDictType model_coeffs() const override;

private:
    //! Kosovic coefficient (default value set for ABL simulations)
    // Ref: Mirocha et. al "Implementation of a Nonlinear Subfilter Turbulence
    //  Stress Model for Large-Eddy Simulation in the Advanced Research WRF
    //  Model"
    // , MWR 2012.
    amrex::Real m_Cb{0.36_rt};
    amrex::Real m_Cs{0.135_rt};
    amrex::Real m_C1{2.1_rt};
    amrex::Real m_C2{2.1_rt};
    amrex::Real m_Sk{0.5_rt};
    bool m_surfaceRANS{false};
    amrex::Real m_surfaceFactor{0.0_rt};
    amrex::Real m_switchLoc{24.0_rt};
    amrex::Real m_surfaceRANSExp{2.0_rt};
    amrex::Real m_LESTurnOff{1.0e15_rt};
    bool m_writeTerms{false};
    const Field& m_vel;
    const Field& m_rho;
    const Field& m_theta;
    Field& m_Nij;
    Field& m_divNij;
    std::string m_wall_het_model{"none"};
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
    amrex::Real m_kappa{0.41};
    amrex::Real m_gamma_m{5.0};
    amrex::Real m_beta_m{16.0};
    amrex::Real m_surface_roughness_z0{0.1};
    amrex::Real m_ref_temp{300.0};
};

} // namespace amr_wind::turbulence

#endif /* SMAGORINSKY_H */
