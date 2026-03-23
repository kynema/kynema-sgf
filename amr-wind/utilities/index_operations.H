#ifndef INDEX_OPERATIONS_H
#define INDEX_OPERATIONS_H
#include <algorithm>

#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "AMReX_Box.H"
#include "AMReX_RealBox.H"
#include "AMReX_RealVect.H"
#include "AMReX_Geometry.H"

using namespace amrex::literals;

namespace amr_wind::utils {

//! Return closest index (from lower) of value in vector
AMREX_FORCE_INLINE int closest_index(
    const amrex::Vector<amrex::Real>& vec,
    const amrex::Real value,
    const amrex::Real tol = 0.0_rt)
{
    auto it = std::ranges::upper_bound(vec, value);
    if (it == vec.end()) {
        // Try again with tolerance added
        it = std::ranges::upper_bound(vec, value - tol);
    }
    AMREX_ALWAYS_ASSERT(it != vec.end());

    const int idx = static_cast<int>(std::distance(vec.begin(), it));
    return amrex::max(idx - 1, 0);
}

//! Return indices perpendicular to normal
template <typename T = amrex::GpuArray<int, 2>>
AMREX_FORCE_INLINE T perpendicular_idx(const int normal)
{
    switch (normal) {
    case 0:
        return T{1, 2};
    case 1:
        return T{0, 2};
    case 2:
        return T{0, 1};
    default:
        amrex::Abort("Invalid normal value to determine perpendicular indices");
    }
    return T{-1, -1};
}

/** Get the intersection with a boundary box while considering if on face or
 * cell. Intended for auxiliary boundary fill calls.
 *
 *  \param shift_to_cc output index vector to shift from output box to
 *  a cell-centered representation.
 *  \param grown_interior_box box grown from domain interior to overlap with
 * boundary
 *  \param domain_boundary_box box representing domain boundary containing data
 *  \param ori orientation of current boundary condition being evaluated
 *  \return The Box instance representing the intersection of the two inputs
 * with shifts added to consider the location of the field associated with the
 * grown_interior_box
 */
AMREX_INLINE amrex::Box face_aware_boundary_box_intersection(
    amrex::IntVect& shift_to_cc,
    amrex::Box grown_interior_box,
    const amrex::Box& domain_boundary_box,
    const amrex::Orientation& ori)
{
    shift_to_cc = {0, 0, 0};
    const auto& field_location_vector = grown_interior_box.type();
    if (!grown_interior_box.cellCentered()) {
        grown_interior_box.enclosedCells();
    }
    auto bx = grown_interior_box & domain_boundary_box;
    if (bx.isEmpty()) {
        return bx;
    }

    if (ori.isHigh() && field_location_vector[ori.coordDir()] == 1) {
        bx.shift(field_location_vector);
        shift_to_cc = -field_location_vector;
    }
    return bx;
}

// This version omits the shift to cell-centered data, often not needed
AMREX_INLINE amrex::Box face_aware_boundary_box_intersection(
    amrex::Box grown_interior_box,
    const amrex::Box& domain_boundary_box,
    const amrex::Orientation& ori)
{
    amrex::IntVect shift_to_cc = {0, 0, 0};
    return face_aware_boundary_box_intersection(
        shift_to_cc, grown_interior_box, domain_boundary_box, ori);
}

//! Check if a point is inside a box
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool contains(
    const amrex::Box& box,
    const amrex::RealVect& pos,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& prob_lo,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxinv)
{
    const amrex::IntVect iv(AMREX_D_DECL(
        static_cast<int>(amrex::Math::floor((pos[0] - prob_lo[0]) * dxinv[0])),
        static_cast<int>(amrex::Math::floor((pos[1] - prob_lo[1]) * dxinv[1])),
        static_cast<int>(
            amrex::Math::floor((pos[2] - prob_lo[2]) * dxinv[2]))));
    return box.contains(iv);
}

/** Convert a bounding box into amrex::Box index space at a given level
 *
 *  \param rbx Bounding box as defined in global domain coordinates
 *  \param geom AMReX geometry information for a given level
 *  \return The Box instance that defines the index space equivalent to bounding
 * box
 */
amrex::Box
realbox_to_box(const amrex::RealBox& rbx, const amrex::Geometry& geom);

} // namespace amr_wind::utils

#endif /* INDEX_OPERATIONS_H */
