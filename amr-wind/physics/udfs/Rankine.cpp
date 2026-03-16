#include "amr-wind/physics/udfs/Rankine.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/vs/vector.H"
#include "amr-wind/equation_systems/icns/icns.H"

#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::udf {

Rankine::Rankine(const Field& fld)
{
    AMREX_ALWAYS_ASSERT(fld.name() == pde::ICNS::var_name());
    AMREX_ALWAYS_ASSERT(fld.num_comp() == AMREX_SPACEDIM);

    const int ncomp = fld.num_comp();

    {
        amrex::ParmParse pp("incflo");
        amrex::Vector<amrex::Real> vel(ncomp, 0.0_rt);
        pp.getarr("velocity", vel);
        AMREX_ALWAYS_ASSERT(vel.size() == ncomp);
        for (int i = 0; i < ncomp; ++i) {
            m_op.vel_ref[i] = vel[i];
        }
    }
    {
        amrex::ParmParse pp("Rankine");
        pp.query("Umax", m_op.Umax);
        pp.query("Rmax", m_op.Rmax);
        amrex::Vector<amrex::Real> start_location;
        pp.queryarr("start_location", start_location);
        for (int i = 0; i < start_location.size(); ++i) {
            m_op.start_location[i] = start_location[i];
        }
    }
}

} // namespace amr_wind::udf
