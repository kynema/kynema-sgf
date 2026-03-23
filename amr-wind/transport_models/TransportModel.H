#ifndef TRANSPORTMODEL_H
#define TRANSPORTMODEL_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/CFDSim.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/utilities/constants.H"
#include "amr-wind/core/field_ops.H"

namespace amr_wind::transport {

/**
 *  \defgroup transport Transport models
 *
 *  AMR-Wind representation of fluid transport models.
 *
 *  For most land-based wind-energy applications, a constant transport
 *  (ConstTransport) model is sufficient with uniform properties. However, for
 *  offshore applications with multiphase models, the properties can change at
 *  different positions in the domain. Transport models provide a way to
 *  separate out this behavior from the base turbulence model implementation.
 *  The turbulence models are templated on transport model to provide different
 *  behavior depending on fluid transport model chosen by the user. For more
 *  details see amr_wind::turbulence::TurbModel implementation.
 *
 */

/** Abstract representation of a transport model
 *  \ingroup transport
 */
class TransportModel : public Factory<TransportModel, const CFDSim&>
{
public:
    static std::string base_identifier() { return "TransportModel"; }

    ~TransportModel() override = default;

    //! Dynamic laminar viscosity (kg/m/s)
    virtual std::unique_ptr<ScratchField> mu() = 0;

    //! Thermal diffusivity
    virtual std::unique_ptr<ScratchField> alpha() = 0;

    //! Scalar diffusivity based on Schmidt number
    virtual std::unique_ptr<ScratchField>
    scalar_diffusivity(const std::string& scalar_name) = 0;

    //! Thermal expansion coefficient
    [[nodiscard]] virtual std::unique_ptr<ScratchField> beta() const = 0;

    //! Thermal expansion coefficient
    virtual void beta_impl(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const amrex::Array4<amrex::Real>& beta) const = 0;

    //! Reference temperature
    [[nodiscard]] virtual amrex::Real reference_temperature() const = 0;

    //! Reference temperature
    virtual void ref_theta_impl(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const amrex::Array4<amrex::Real>& ref_theta) const = 0;

    //! Reference temperature
    [[nodiscard]] virtual std::unique_ptr<ScratchField> ref_theta() const = 0;

protected:
    const int m_ngrow{1};
};
} // namespace amr_wind::transport

#endif /* TRANSPORTMODEL_H */
