#ifndef PDEFIELDS_H
#define PDEFIELDS_H

#include <iosfwd>

namespace amr_wind {

class FieldRepo;
class Field;

namespace pde {

/** Collection of fields related to a PDE
 *  \ingroup eqsys
 *
 *  Every PDE requires a common set of fields, e.g., the transported variable
 *  (scalar or vector), advection, diffusion, source, and effective viscosity
 *  fields. PDEFields is a struct that holds these fields and provides a common
 *  way to access these fields regardless of the PDE that is being solved.
 *
 *  The fields are all named in a consistent manner based on the variable name
 *  as defined by the PDE trait (e.g., ICNS, Temperature, etc.). For example, if
 *  the user is solving the TKE equation, then the field variable is `tke` and
 *  the source term is `tke_src_term`, and so on.
 */
struct PDEFields
{
    PDEFields(FieldRepo& repo_in, const std::string& var_name);

    //! Reference to the field repository instance
    FieldRepo& repo;

    //! Solution variable (e.g., velocity, temperature)
    Field& field;
    //! Effective visocity field (e.g., velocity_mueff)
    Field& mueff;

    //! Source term for this PDE
    Field& src_term;
    //! Diffusion term for this PDE
    Field& diff_term;
    //! Convective term for this PDE
    Field& conv_term;
};

} // namespace pde
} // namespace amr_wind

#endif /* PDEFIELDS_H */
