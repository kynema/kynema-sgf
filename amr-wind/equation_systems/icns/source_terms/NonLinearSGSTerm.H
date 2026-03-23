#ifndef NONLINEARSGSTERM_H
#define NONLINEARSGSTERM_H

#include "amr-wind/equation_systems/icns/MomentumSource.H"
#include "amr-wind/core/SimTime.H"

namespace amr_wind::pde::icns {

/** Adds the non-linear terms from the Kosovic Model to Momentum Equation..
 *
 *  \ingroup icns_src NonLinearSGS
 *
 *  \sa NonLinearSGS
 */
class NonLinearSGSTerm : public MomentumSource::Register<NonLinearSGSTerm>
{
public:
    static std::string identifier() { return "NonLinearSGSTerm"; }

    explicit NonLinearSGSTerm(const CFDSim& sim);

    ~NonLinearSGSTerm() override;

    void operator()(
        int lev,
        const amrex::MFIter& mfi,
        const amrex::Box& bx,
        FieldState fstate,
        const amrex::Array4<amrex::Real>& src_term) const override;

private:
    const Field& m_divNij;
};

} // namespace amr_wind::pde::icns

#endif /* NONLINEARSGSTERM_H */
