#include "amr-wind/utilities/index_operations.H"

namespace amr_wind::utils {
amrex::Box
realbox_to_box(const amrex::RealBox& rbx, const amrex::Geometry& geom)
{
    const auto* problo = geom.ProbLo();
    const auto* probhi = geom.ProbHi();
    const auto* dxi = geom.InvCellSize();

    amrex::IntVect lo, hi;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        amrex::Real bbox_min = amrex::max(rbx.lo()[i], problo[i]);
        amrex::Real bbox_max = amrex::min(rbx.hi()[i], probhi[i]);

        amrex::Real rlo = std::floor((bbox_min - problo[i]) * dxi[i]);
        amrex::Real rhi = std::ceil((bbox_max - problo[i]) * dxi[i]);

        lo[i] = static_cast<int>(rlo);
        hi[i] = static_cast<int>(rhi);
    }

    return amrex::Box{lo, hi};
}

} // namespace amr_wind::utils
