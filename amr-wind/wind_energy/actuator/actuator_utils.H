#ifndef ACUTATOR_UTILS_H
#define ACUTATOR_UTILS_H

#include <set>
#include <cmath>
#include "amr-wind/core/vs/vector_space.H"
#include "amr-wind/utilities/index_operations.H"
#include "AMReX_AmrCore.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator {

struct ActInfo;

namespace utils {

/** Return a set of process IDs (MPI ranks) that contain AMR boxes that interact
 *  with a given actuator body.
 *
 *  The region of influence of the actuator body is determined by checking for
 *  intersections with a bounding box ``rbox``.
 *
 *  \param mesh AMReX mesh instance
 *  \param rbx The bounding box that defines the region of influence of a
 * turbine
 */
std::set<int> determine_influenced_procs(
    const amrex::AmrCore& mesh, const amrex::RealBox& rbx);

void determine_root_proc(
    ActInfo& /*info*/, amrex::Vector<int>& /*act_proc_count*/);

/** Return the Gaussian smearing factor in 3D
 *
 *  \param dist Distance vector of the cell center from the actuator node in
 *  local frame of reference.
 *
 *  \param eps Three-dimensional Gaussian scaling factor
 *
 *  \return Gaussian smearing factor in 3D
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
gaussian3d(const vs::Vector& dist, const vs::Vector& eps)
{
    const vs::Vector rr{
        dist.x() / eps.x(), dist.y() / eps.y(), dist.z() / eps.z()};
    const amrex::Real rr_sqr = vs::mag_sqr(rr);

    if (rr_sqr < 16.0_rt) {
        constexpr amrex::Real fac = 0.17958712212516656_rt;
        const amrex::Real eps_fac = eps.x() * eps.y() * eps.z();
        return (fac / eps_fac) *
               std::exp(
                   (-rr.x() * rr.x()) - (rr.y() * rr.y()) - (rr.z() * rr.z()));
    }

    return 0.0_rt;
}
/** Return the Gaussian smearing factor in 1D
 *
 *  \param dist Distance vector of the cell center from the actuator node in
 *  local frame of reference.
 *
 *  \param eps One-dimensional Gaussian scaling factor
 *
 *  \return Gaussian smearing factor in 1D
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
gaussian1d(const amrex::Real& dist, const amrex::Real& eps)
{
    const amrex::Real fac = std::numbers::inv_sqrtpi_v<amrex::Real>;
    if (std::abs(dist / eps) >= 16.0_rt) {
        return 0.0_rt;
    }
    return fac / eps * std::exp(-(dist * dist) / (eps * eps));
}

/** Compute distance between two points
 *
 * \param origin Origin for the cylindrical coordinate system
 * \param normal Normal  relative to origin (defines cylinder orientation)
 * \param point1 location of the first point in Cartesian coordinates
 * \param point2 location of second point in Cartesian coordinates
 * @return distance between points in cylindrical coordinates
 * (angle is in radians)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE vs::Vector delta_pnts_cyl(
    const vs::Vector& origin,
    const vs::Vector& normal,
    const vs::Vector& point1,
    const vs::Vector& point2)
{
    const amrex::Real eps = vs::DTraits<amrex::Real>::eps();
    const auto d1 = point1 - origin;
    const auto d2 = point2 - origin;

    // cppcheck-suppress duplicateExpression
    const amrex::Real norm_dist1 = (d1 & normal) / (normal & normal);
    // cppcheck-suppress duplicateExpression
    const amrex::Real norm_dist2 = (d2 & normal) / (normal & normal);

    const auto proj1 = normal * norm_dist1;
    const auto proj2 = normal * norm_dist2;

    // compute vector on disk via vector rejection
    const auto v1 = d1 - proj1;
    const auto v2 = d2 - proj2;

    // can do magnitude since we've removed the normal component
    const amrex::Real r1 = vs::mag(v1);
    const amrex::Real r2 = vs::mag(v2);

    const amrex::Real theta =
        std::acos(amrex::min<amrex::Real>((v1 & v2) / (r1 * r2 + eps), 1.0_rt));

    return {std::abs(r1 - r2), theta, norm_dist1 - norm_dist2};
}

/** Compute weighting for a normalized linear interpolant
 *
 *
 * \param distance distance between sample location and the reference points
 * \param dX spacing between reference points
 * @return normalized integration weight at the sample location
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
linear_basis_1d(const amrex::Real distance, const amrex::Real dX)
{
    return amrex::max<amrex::Real>(0.0_rt, 1.0_rt - (distance / dX)) / dX;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE vs::Vector compute_tangential_vector(
    const vs::Vector& center, const vs::Vector& normal, const vs::Vector& point)
{
    const auto blade_axis = (point - center).normalize();
    return vs::quaternion(normal, 90.0_rt) & blade_axis;
}

} // namespace utils
} // namespace amr_wind::actuator

#endif /* ACUTATOR_UTILS_H */
