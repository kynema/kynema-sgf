#include "amr-wind/physics/udfs/LinearProfile.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::udf {

LinearProfile::LinearProfile(const Field& fld)
{
    amrex::ParmParse pp("LinearProfile." + fld.name());

    const int ncomp = fld.num_comp();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        (ncomp <= AMREX_SPACEDIM),
        "LinearProfile requires field with 3 or fewer components");
    pp.query("direction", m_op.idir);
    const auto geom = fld.repo().mesh().Geom(0);
    m_op.zmin = geom.ProbLo(m_op.idir);
    m_op.zmax = geom.ProbHi(m_op.idir);
    pp.query("start", m_op.zmin);
    pp.query("stop", m_op.zmax);
    amrex::Vector<amrex::Real> start_val, end_val;
    pp.getarr("start_val", start_val);
    pp.getarr("stop_val", end_val);

    AMREX_ALWAYS_ASSERT(start_val.size() == ncomp);
    AMREX_ALWAYS_ASSERT(end_val.size() == ncomp);

    for (int i = 0; i < ncomp; ++i) {
        m_op.fld_min[i] = start_val[i];
        m_op.fld_max[i] = end_val[i];
    }
}

} // namespace amr_wind::udf
