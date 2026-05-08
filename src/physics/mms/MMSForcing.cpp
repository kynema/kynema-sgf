#include "src/physics/mms/MMSForcing.H"
#include "src/CFDSim.H"
#include "masa.h"

namespace kynema_sgf::pde::icns::mms {

/** MMS forcing term
 */
MMSForcing::MMSForcing(const CFDSim& sim)
    : m_mms_vel_source(sim.repo().get_field("mms_vel_source"))
{
    static_assert(AMREX_SPACEDIM == 3, "MMS implementation requires 3D domain");
}

void MMSForcing::operator()(
    const int lev, const FieldState /*fstate*/, amrex::MultiFab& src_term) const
{
    auto const& src_arrs = src_term.arrays();
    auto const& mms_arrs = m_mms_vel_source(lev).const_arrays();

    amrex::ParallelFor(
        src_term, amrex::IntVect(0), AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) {
            src_arrs[nbx](i, j, k, n) += mms_arrs[nbx](i, j, k, n);
        });
}
} // namespace kynema_sgf::pde::icns::mms
