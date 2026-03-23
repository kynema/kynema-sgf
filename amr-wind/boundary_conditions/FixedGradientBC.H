#ifndef FIXEDGRADIENTBC_H
#define FIXEDGRADIENTBC_H

#include "amr-wind/core/FieldBCOps.H"
#include "amr-wind/core/FieldRepo.H"
#include "AMReX_Orientation.H"

namespace amr_wind {

/** Applies a non-zero gradient at the domain boundary
 *  \ingroup field_bc
 *
 *  Used with the inhomogeneous Neumann linear solver BC to apply a
 *  user-specified fixed gradient at the domain boundary.
 */
class FixedGradientBC : public FieldBCIface
{
public:
    FixedGradientBC(Field& field, amrex::Orientation ori);

    void operator()(Field& field, FieldState /*rho_state*/) override;

private:
    Field& m_field;

    amrex::Orientation m_ori;
};

} // namespace amr_wind

#endif /* FIXEDGRADIENTBC_H */
