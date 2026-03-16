#include "amr-wind/physics/udfs/CustomScalar.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/vs/vector.H"
#include "amr-wind/equation_systems/icns/icns.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::udf {

CustomScalar::CustomScalar(const Field& fld)
{
    // This is a where the user can set some user defined variables
    // This capability can be activated with the following in the input file:
    // xlo.type = "mass_inflow"
    // xlo.temperature.inflow_type = CustomScalar
    // CustomScalar.foo = 1.0_rt

    amrex::ParmParse pp("CustomScalar");
    pp.query("foo", m_op.foo);
    const int ncomp = fld.num_comp();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        (ncomp == 1), "CustomScalar requires field with 1 component");

    amrex::Abort(
        "Please define the body of this function and the corresponding struct "
        "in the header file before using it. Then remove this message");
}

} // namespace amr_wind::udf
