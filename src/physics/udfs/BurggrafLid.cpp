#include "src/physics/udfs/BurggrafLid.H"
#include "src/core/Field.H"
#include "src/core/FieldRepo.H"
#include "src/core/vs/vector.H"
#include "src/equation_systems/icns/icns.H"

#include "AMReX_ParmParse.H"

namespace kynema_sgf::udf {

BurggrafLid::BurggrafLid(const Field& fld)
{
    AMREX_ALWAYS_ASSERT(fld.name() == pde::ICNS::var_name());
    AMREX_ALWAYS_ASSERT(fld.num_comp() == AMREX_SPACEDIM);
}

} // namespace kynema_sgf::udf
