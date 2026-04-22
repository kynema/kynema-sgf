#ifndef MASSINFLOWOUTFLOWBC_H
#define MASSINFLOWOUTFLOWBC_H

#include "src/core/FieldBCOps.H"
#include "src/core/FieldRepo.H"
#include "AMReX_Orientation.H"

namespace kynema_sgf {

/** Custom Neumann fills for the inflow-outflow BC
 *  \ingroup field_bc
 *
 *  Used to fill the outflow boundary cells for mass-inflow-outflow BC
 */
class MassInflowOutflowBC : public FieldBCIface
{
public:
    MassInflowOutflowBC(Field& field, amrex::Orientation ori);

    void operator()(Field& field, FieldState /*rho_state*/) override;

private:
    Field& m_field;

    amrex::Orientation m_ori;
};

} // namespace kynema_sgf

#endif /* MASSINFLOWOUTFLOWBC_H */
