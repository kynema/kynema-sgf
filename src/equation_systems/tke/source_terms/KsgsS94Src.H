#ifndef KSGSS94SRC_H
#define KSGSS94SRC_H

#include "src/equation_systems/tke/TKESource.H"

namespace kynema_sgf {
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

    void operator()(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override;
};
} // namespace tke
} // namespace pde
} // namespace kynema_sgf

#endif /* KSGSS94SRC_H */
