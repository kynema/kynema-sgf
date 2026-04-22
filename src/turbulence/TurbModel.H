#ifndef TURBMODEL_H
#define TURBMODEL_H

#include "amr-wind/turbulence/TurbulenceModel.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::turbulence {

/** Turbulence model with a specific transport model
 *  \ingroup turb_model
 *
 *  This templated class allows turbulence model to be specialised based on \ref
 *  transport models. All turbulence model implementations inherit from this
 * class.
 */
template <typename Transport>
class TurbModel : public TurbulenceModel
{
public:
    using TransportType = Transport;

    explicit TurbModel(CFDSim& sim) : m_sim(sim), m_transport(sim) {}

    void register_mueff_field(Field& mueff) override { m_mueff = &mueff; }

    void register_alphaeff_field(Field& alphaeff) override
    {
        m_alphaeff = &alphaeff;
    }

    //! Return the dynamic viscosity field
    std::unique_ptr<ScratchField> mu() override { return m_transport.mu(); }

    //! Return the thermal diffusivity field
    std::unique_ptr<ScratchField> alpha() override
    {
        return m_transport.alpha();
    }

    //! Return the scalar diffusivity field
    std::unique_ptr<ScratchField>
    scalar_diffusivity(const std::string& name) override
    {
        return m_transport.scalar_diffusivity(name);
    }

    Field& mueff() override
    {
        AMREX_ASSERT(m_mueff != nullptr);
        return *m_mueff;
    }

    Field& alphaeff() override
    {
        AMREX_ASSERT(m_alphaeff != nullptr);
        return *m_alphaeff;
    }

protected:
    CFDSim& m_sim;

    //! Transport properties instance
    Transport m_transport;

    //! Reference to the effective viscosity field
    Field* m_mueff{nullptr};

    //! Reference to the thermal diffusivity field
    Field* m_alphaeff{nullptr};
};

} // namespace amr_wind::turbulence

#endif /* TURBMODEL_H */
