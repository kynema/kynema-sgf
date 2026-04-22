#ifndef LAMINARMODEL_H
#define LAMINARMODEL_H

#include "amr-wind/turbulence/TurbModel.H"

namespace amr_wind::turbulence {

/** Laminar model
 *  \ingroup turb_model
 *
 *  Implementation of laminar flow model
 */
template <typename Transport>
class Laminar : public TurbModel<Transport>
{
public:
    static std::string identifier()
    {
        return "Laminar-" + Transport::identifier();
    }

    explicit Laminar(CFDSim& sim) : TurbModel<Transport>(sim) {}

    ~Laminar() override = default;

    //! String representation of this model for I/O
    [[nodiscard]] std::string model_name() const override { return "laminar"; }

    //! Update the effective/turbulent viscosity field
    void update_turbulent_viscosity(
        FieldState fstate, DiffusionType diff_type /*unused*/) override;

    //! No post advance work for this model
    void post_advance_work() override {}

    //! No model specific coefficient for this model
    void parse_model_coeffs() override {}

    //! Return the turbulent viscosity field (just the same as the effective
    //! field)
    Field& mu_turb() override { return this->mueff(); }

    //! Return the thermal diffusivity field (just the same as effective field)
    Field& alpha_turb() override { return this->alphaeff(); }

    //! Indicate that this model is not a turbulent model type
    [[nodiscard]] bool is_turbulent() const override { return false; }

    //! Interface to update effective viscosity (mu_eff = mu + mu_t)
    void update_mueff(Field& mueff) override;

    //! Interface to update effective thermal diffusivity
    void update_alphaeff(Field& alphaeff) override;

    //! Interface to update scalar diffusivity based on Schmidt number
    void update_scalar_diff(Field& deff, const std::string& name) override;

    //! Return model coefficients dictionary
    [[nodiscard]] TurbulenceModel::CoeffsDictType model_coeffs() const override;
};

} // namespace amr_wind::turbulence

#endif /* LAMINARMODEL_H */
