#ifndef KWSSTSRC_H
#define KWSSTSRC_H

#include "amr-wind/equation_systems/tke/TKESource.H"

namespace amr_wind::pde::tke {

/** Source terms for TKE equation in K-w-SST turbulence model
 *   \ingroup turb_model
 *
 *  Works for all TKE equations with production and dissipation terms
 */
class KwSSTSrc : public TKESource::Register<KwSSTSrc>
{
public:
    static std::string identifier() { return "KwSSTSrc"; }

    explicit KwSSTSrc(const CFDSim& /*sim*/);

    ~KwSSTSrc() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    Field& m_shear_prod;
    Field& m_diss;
    Field& m_buoy_term;
};

} // namespace amr_wind::pde::tke

#endif /* KWSSTSRC_H */
