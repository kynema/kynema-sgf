#ifndef ONEEQKSGS_H
#define ONEEQKSGS_H

#include <string>
#include "src/core/ScratchField.H"
#include "src/turbulence/TurbModelBase.H"
#include "src/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::turbulence {

/** Base class for 1-Equation subgrid scale TKE turbulence model
 *  \ingroup turb_model
 */
template <typename Transport>
class OneEqKsgs : public TurbModelBase<Transport>
{
public:
    explicit OneEqKsgs(CFDSim& sim);

    ~OneEqKsgs() override;

protected:
    Field& m_vel;

    Field& m_turb_lscale;
    Field& m_shear_prod;
    Field& m_buoy_prod;
    Field& m_dissip;
    Field& m_rho;

    Field* m_tke{nullptr};

    //! Turbulence constant
    amrex::Real m_Ce{0.1_rt};
    amrex::Real m_Ceps{0.93_rt};
};

/** 1-equation subgrid scale TKE turbulence model
 *  \ingroup turb_model we_abl
 *
 *  Implements the model described in Moeng 1984 paper
 *
 *  \sa kynema_sgf::pde::tke::KsgsM84Src
 */
template <typename Transport>
class OneEqKsgsM84 : public OneEqKsgs<Transport>
{
public:
    static std::string identifier()
    {
        return "OneEqKsgsM84-" + Transport::identifier();
    }

    explicit OneEqKsgsM84(CFDSim& sim);

    ~OneEqKsgsM84() override;

    [[nodiscard]] std::string model_name() const override
    {
        return "OneEqKsgsM84";
    }

    //! Update the turbulent viscosity field
    void update_turbulent_viscosity(
        FieldState fstate, DiffusionType /*unused*/) override;

    //! Do any post advance work
    void post_advance_work() override;

    //! Actions to perform post initialization
    void post_init_actions() override;

    //! Actions to perform post regrid
    void post_regrid_actions() override;

    //! Update the effective thermal diffusivity field
    void update_alphaeff(Field& alphaeff) override;

    //! Update the effective scalar diffusivity field
    void update_scalar_diff(Field& deff, const std::string& name) override;

    //! Parse turbulence model coefficients
    void parse_model_coeffs() override;

    //! Return turbulence model coefficients
    [[nodiscard]] TurbulenceModel::CoeffsDictType model_coeffs() const override;

    //! Update the near-wall hybrid length scale field (surfaceRANS option)
    void update_hybrid_length_scale();

private:
    Field& m_temperature;

    //! Gravity vector (m/s^2)
    amrex::Vector<amrex::Real> m_gravity{0.0_rt, 0.0_rt, -9.81_rt};

    //! Hybrid RANS-LES with Nalu-wind
    bool m_hybrid_rl{false};
    Field* m_sdr;

    std::unique_ptr<ScratchField> m_gradT;

    /** Near-wall hybrid length scale (OneEqKsgsM84.surfaceRANS)
     *
     *  Blends the filter width with the wall-distance length scale that
     *  recovers the log law, l^2 = (1-f)^n ds^2 + f^n l_rans^2 with
     *  f = exp(-z / switchLoc), so that the length scale near the ground does
     *  not depend on the refinement level. Same blending and inputs as the
     *  Kosovic model (surfaceRANS, switchLoc, surfaceRANSExp).
     */
    bool m_surface_rans{false};
    amrex::Real m_switch_loc{24.0_rt};
    amrex::Real m_surface_rans_exp{2.0_rt};
    amrex::Real m_kappa{0.41_rt};
    amrex::Real m_gamma_m{5.0_rt};
    amrex::Real m_beta_m{16.0_rt};
    //! Obukhov length: ABL.monin_obukhov_length when given, else the value
    //! of the ABL wall function at the current step
    amrex::Real m_monin_obukhov_length{constants::LARGE_NUM};
    bool m_fixed_mol{false};
    Field* m_hybrid_lscale{nullptr};
};

/** 1-equation subgrid scale TKE turbulence model
 *  \ingroup turb_model we_abl
 *
 *  Implements the model described in Sullivan 1994 paper
 *
 *  \sa kynema_sgf::pde::tke::KsgsS94Src
 */
template <typename Transport>
class OneEqKsgsS94 : public OneEqKsgs<Transport>
{
public:
    static std::string identifier()
    {
        return "OneEqKsgsS94-" + Transport::identifier();
    }

    explicit OneEqKsgsS94(CFDSim& sim);

    ~OneEqKsgsS94() override;

    [[nodiscard]] std::string model_name() const override
    {
        return "OneEqKsgsS94";
    }

    //! Update the turbulent viscosity field
    void update_turbulent_viscosity(
        FieldState fstate, DiffusionType diff_type /*unused*/) override;

    //! No post advance work for this model
    void post_advance_work() override {}

    //! Update the effective scalar diffusivity field
    void update_scalar_diff(Field& deff, const std::string& name) override;

    //! Parse turbulence model coefficients
    void parse_model_coeffs() override;

    //! Return turbulence model coefficients
    [[nodiscard]] TurbulenceModel::CoeffsDictType model_coeffs() const override;
};

} // namespace kynema_sgf::turbulence

#endif /* ONEEQKSGS_H */
