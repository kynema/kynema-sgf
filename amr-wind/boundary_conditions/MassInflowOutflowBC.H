#ifndef MASSINFLOWOUTFLOWBC_H
#define MASSINFLOWOUTFLOWBC_H

#include "amr-wind/core/FieldBCOps.H"
#include "amr-wind/core/FieldRepo.H"
#include "AMReX_Orientation.H"

namespace amr_wind {

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

} // namespace amr_wind

#endif /* MASSINFLOWOUTFLOWBC_H */
