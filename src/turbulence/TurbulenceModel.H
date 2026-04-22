#ifndef TURBULENCEMODEL_H
#define TURBULENCEMODEL_H

#include <unordered_map>

#include "amr-wind/core/Factory.H"
#include "amr-wind/core/FieldDescTypes.H"
#include "amr-wind/incflo.H"

namespace amr_wind {

class CFDSim;
class ScratchField;
class Field;

namespace turbulence {

/**
 *  \defgroup turb_model Turbulence models
 *
 *  Turbulence models available in AMR-Wind
 */

/** Abstract representation of a turbulence model
 *  \ingroup turb_model
 *
 *  This class is an abstract representation of a turbulence model used in
 *  AMR-Wind. Subclasses must implement the methods defined in this class to
 *  interact with \ref eqsys.
 */
class TurbulenceModel : public Factory<TurbulenceModel, CFDSim&>
{
public:
    //! Default identifier used by factor methods
    using CoeffsDictType = std::unordered_map<std::string, amrex::Real>;

    static std::string base_identifier() { return "TurbulenceModel"; }

    ~TurbulenceModel() override = default;

    /** Update the turbulent viscosity field
     *
     *  This method is called at the beginning of each predictor/corrector step
     *  to update the turbulent viscosity field. The exact effective viscosity
     *  for each equation system is handled by the concrete implementation of
     *  PDESystem.
     *
     *  \param fstate State used for updates to differentiate logic in predictor
     *  and corrector steps.
     *  \param diff_type diffusion solve type
     */
    virtual void
    update_turbulent_viscosity(FieldState fstate, DiffusionType diff_type) = 0;

    //! Do any post advance actions for the turbulence model
    virtual void post_advance_work() = 0;

    //! Register the effective viscosity (for momentum field)
    virtual void register_mueff_field(Field& mueff) = 0;

    //! Register the thermal diffusivity field
    virtual void register_alphaeff_field(Field& alphaeff) = 0;

    //! String representation of this turbulence model for I/O
    [[nodiscard]] virtual std::string model_name() const = 0;

    //! Return the dynamic viscosity (laminar) field
    virtual std::unique_ptr<ScratchField> mu() = 0;

    //! Return the thermal diffusivity (laminar) field for enthalpy/temperature
    virtual std::unique_ptr<ScratchField> alpha() = 0;

    //! Return the scalar diffusivity field
    virtual std::unique_ptr<ScratchField>
    scalar_diffusivity(const std::string& name) = 0;

    //! Return the turbulent dynamic viscosity field
    virtual Field& mu_turb() = 0;

    //! Return the turbulent dynamic viscosity field
    virtual Field& alpha_turb() = 0;

    //! Return the effective dynamic viscosity field
    virtual Field& mueff() = 0;

    //! Return the effective thermal diffusivity field
    virtual Field& alphaeff() = 0;

    //! Flag indicating whether the model is turbulent
    [[nodiscard]] virtual bool is_turbulent() const { return true; }

    //! Interface to update effective viscosity
    //!
    //! \f$\mu_\mathrm{eff} = \mu + \mu_t\f$
    virtual void update_mueff(Field& mueff) = 0;

    //! Interface to update effective thermal diffusivity
    virtual void update_alphaeff(Field& alphaeff) = 0;

    //! Interface to update scalar diffusivity based on Schmidt number
    virtual void update_scalar_diff(Field& deff, const std::string& name) = 0;

    //! Parse turbulence model coefficients
    virtual void parse_model_coeffs() = 0;

    //! Return model coefficients dictionary
    [[nodiscard]] virtual CoeffsDictType model_coeffs() const = 0;
};

} // namespace turbulence
} // namespace amr_wind

#endif /* TURBULENCEMODEL_H */
