#include "src/physics/udfs/CustomVelocity.H"
#include "src/core/Field.H"
#include "src/core/FieldRepo.H"
#include "src/core/vs/vector.H"
#include "src/equation_systems/icns/icns.H"

#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::udf {

CustomVelocity::CustomVelocity(const Field& fld)
{
    // This is a where the user can set some user defined variables
    // This capability can be activated with the following in the input file:
    // xlo.type = "mass_inflow"
    // xlo.velocity.inflow_type = CustomVelocity
    // CustomVelocity.foo = 1.0

    const int ncomp = fld.num_comp();
    amrex::ParmParse pp("CustomVelocity");
    pp.query("foo", m_op.foo);
    amrex::Vector<amrex::Real> vel(ncomp, 0.0_rt);
    pp.getarr("velocity", vel);
    AMREX_ALWAYS_ASSERT(vel.size() == ncomp);
    for (int i = 0; i < ncomp; ++i) {
        m_op.bar[i] = vel[i];
    }
    amrex::Abort(
        "Please define the body of this function and the corresponding struct "
        "in the header file before using it. Then remove this message");
}

} // namespace kynema_sgf::udf
