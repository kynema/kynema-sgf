#include "amr-wind/physics/mms/MMSForcing.H"
#include "amr-wind/CFDSim.H"
#include "masa.h"

namespace amr_wind::pde::icns::mms {

/** MMS forcing term
 */
MMSForcing::MMSForcing(const CFDSim& sim)
    : m_mms_vel_source(sim.repo().get_field("mms_vel_source"))
{
    static_assert(AMREX_SPACEDIM == 3, "MMS implementation requires 3D domain");
}

void MMSForcing::operator()(
    const int lev,
    const amrex::MFIter& mfi,
    const amrex::Box& bx,
    const FieldState /*fstate*/,
    const amrex::Array4<amrex::Real>& src_term) const
{
    const auto& mms_src = m_mms_vel_source(lev);
    const auto& mms_src_arr = mms_src.const_array(mfi);

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        src_term(i, j, k, 0) += mms_src_arr(i, j, k, 0);
        src_term(i, j, k, 1) += mms_src_arr(i, j, k, 1);
        src_term(i, j, k, 2) += mms_src_arr(i, j, k, 2);
    });
}
} // namespace amr_wind::pde::icns::mms
