#include "amr-wind/physics/udfs/PowerLawProfile.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/vs/vector.H"
#include "amr-wind/equation_systems/icns/icns.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::udf {

PowerLawProfile::PowerLawProfile(const Field& fld)
{
    AMREX_ALWAYS_ASSERT(fld.name() == pde::ICNS::var_name());
    AMREX_ALWAYS_ASSERT(fld.num_comp() == AMREX_SPACEDIM);

    amrex::ParmParse pp("PowerLawProfile");
    pp.query("direction", m_op.idir);
    amrex::Vector<amrex::Real> vel;
    pp.get("zref", m_op.zref);
    pp.get("shear_exponent", m_op.shear_exp);
    pp.getarr("uref", vel);
    pp.query("zoffset", m_op.zoffset);
    pp.query("umin", m_op.umin);
    pp.query("umax", m_op.umax);

    AMREX_ALWAYS_ASSERT(vel.size() == AMREX_SPACEDIM);
    m_op.uref = vs::mag(vs::Vector{vel[0], vel[1], vel[2]});
    m_op.umin /= m_op.uref;
    m_op.umax /= m_op.uref;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        m_op.uvec[i] = vel[i];
    }
}

} // namespace amr_wind::udf
