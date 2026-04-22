#ifndef SDRSRC_H
#define SDRSRC_H

#include "amr-wind/equation_systems/sdr/SDRSource.H"

namespace amr_wind::pde::tke {

/** Source terms for SDR equation
 *   \ingroup turb_model
 *
 *  Works for all SDR equations with source term
 */
class SDRSrc : public SDRSource::Register<SDRSrc>
{
public:
    static std::string identifier() { return "SDRSrc"; }

    explicit SDRSrc(const CFDSim& /*sim*/);

    ~SDRSrc() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    Field& m_sdr_src;
    Field& m_sdr_diss;
};

} // namespace amr_wind::pde::tke

#endif /* SDRSRC_H */
