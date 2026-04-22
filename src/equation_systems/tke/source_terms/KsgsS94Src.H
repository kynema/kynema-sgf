#ifndef KSGSS94SRC_H
#define KSGSS94SRC_H

#include "amr-wind/equation_systems/tke/TKESource.H"

namespace amr_wind {
namespace pde {
namespace tke {

/** One-equation k-sgs source term based on Sullivan 1994 paper
 *  \ingroup tke_src turb_model we_abl
 */
class KsgsS94Src : public TKESource::Register<KsgsS94Src>
{
public:
    static const std::string identifier() { return "KsgsS94Src"; }

    KsgsS94Src(const CFDSim&);

    virtual ~KsgsS94Src();

    virtual void operator()(
        const int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        const FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;
};
} // namespace tke
} // namespace pde
} // namespace amr_wind

#endif /* KSGSS94SRC_H */
