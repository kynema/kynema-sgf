#include "amr-wind/physics/udfs/TwoLayer.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/FieldRepo.H"
#include "amr-wind/core/vs/vector.H"
#include "amr-wind/equation_systems/icns/icns.H"

#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::udf {

TwoLayer::TwoLayer(const Field& fld)
{
    // For a 2-layer flow with a top and a bottom layer
    // divided at a specific z-coordinate (default: 0.5_rt)
    // and an optional initial perturbation (default: 1.0_rt)

    {
        const int ncomp = fld.num_comp();
        amrex::ParmParse pp("TwoLayer");

        amrex::Vector<amrex::Real> top_vel(ncomp, 0.0_rt);
        amrex::Vector<amrex::Real> bottom_vel(ncomp, 0.0_rt);
        pp.getarr("top_vel", top_vel);
        pp.getarr("bottom_vel", bottom_vel);

        AMREX_ALWAYS_ASSERT(top_vel.size() == ncomp);
        AMREX_ALWAYS_ASSERT(bottom_vel.size() == ncomp);
        for (int i = 0; i < ncomp; ++i) {
            m_op.top_vel[i] = top_vel[i];
            m_op.bottom_vel[i] = bottom_vel[i];
        }

        pp.query("init_perturb", m_op.init_perturb);
        pp.query("z_partition", m_op.z_part);
    }
}

} // namespace amr_wind::udf
