#ifndef VS_QUATERNION_H
#define VS_QUATERNION_H

#include "src/core/vs/tensor.H"
#include "src/utilities/constants.H"
#include "src/utilities/trig_ops.H"

#include <algorithm>
#include <cmath>
#include <limits>

#include "AMReX.H"

namespace kynema_sgf::vs {

using namespace amrex::literals;

/** Scalar-first unit quaternion used to interpolate and compose rotations. */
struct Quaternion
{
    amrex::Real w{1.0_rt};
    amrex::Real x{0.0_rt};
    amrex::Real y{0.0_rt};
    amrex::Real z{0.0_rt};
};

inline Quaternion normalized(Quaternion q)
{
    const amrex::Real norm =
        std::sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
    if (norm <= std::numeric_limits<amrex::Real>::epsilon()) {
        amrex::Abort("Cannot normalize a zero quaternion");
    }
    q.w /= norm;
    q.x /= norm;
    q.y /= norm;
    q.z /= norm;
    return q;
}

inline amrex::Real dot(const Quaternion& a, const Quaternion& b)
{
    return a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Quaternion conjugate(const Quaternion& q)
{
    return {q.w, -q.x, -q.y, -q.z};
}

//! Rotation composition: apply rhs first, followed by lhs
inline Quaternion operator*(const Quaternion& lhs, const Quaternion& rhs)
{
    return {
        lhs.w * rhs.w - lhs.x * rhs.x - lhs.y * rhs.y - lhs.z * rhs.z,
        lhs.w * rhs.x + lhs.x * rhs.w + lhs.y * rhs.z - lhs.z * rhs.y,
        lhs.w * rhs.y - lhs.x * rhs.z + lhs.y * rhs.w + lhs.z * rhs.x,
        lhs.w * rhs.z + lhs.x * rhs.y - lhs.y * rhs.x + lhs.z * rhs.w};
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor tensor_unit(const Quaternion& q)
{
    return {1.0_rt - 2.0_rt * (q.y * q.y + q.z * q.z),
            2.0_rt * (q.x * q.y - q.z * q.w),
            2.0_rt * (q.x * q.z + q.y * q.w),
            2.0_rt * (q.x * q.y + q.z * q.w),
            1.0_rt - 2.0_rt * (q.x * q.x + q.z * q.z),
            2.0_rt * (q.y * q.z - q.x * q.w),
            2.0_rt * (q.x * q.z - q.y * q.w),
            2.0_rt * (q.y * q.z + q.x * q.w),
            1.0_rt - 2.0_rt * (q.x * q.x + q.y * q.y)};
}

inline Tensor tensor(const Quaternion& input)
{
    return tensor_unit(normalized(input));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Quaternion
from_axis_angle(const Vector& axis, const amrex::Real angle_degrees)
{
    const amrex::Real half_angle = -0.5_rt * utils::radians(angle_degrees);
    const amrex::Real axis_magnitude = mag(axis);
    const amrex::Real scale = std::sin(half_angle) / axis_magnitude;
    return {
        std::cos(half_angle), scale * axis.x(), scale * axis.y(),
        scale * axis.z()};
}

//! Compatibility interface for existing axis-angle tensor rotations
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Tensor
quaternion(const Vector& axis, const amrex::Real angle_degrees)
{
    // Preserve the original operation sequence used by existing actuator and
    // sampling models. Although tensor_unit(from_axis_angle(...)) is
    // algebraically equivalent, it changes floating-point rounding and
    // therefore existing regression results.
    const amrex::Real angle = -1.0_rt * utils::radians(angle_degrees);
    const amrex::Real cval = std::cos(0.5_rt * angle);
    const amrex::Real sval = std::sin(0.5_rt * angle);
    const auto axis_magnitude = mag(axis);
    const amrex::Real q0 = cval;
    const amrex::Real q1 = sval * axis.x() / axis_magnitude;
    const amrex::Real q2 = sval * axis.y() / axis_magnitude;
    const amrex::Real q3 = sval * axis.z() / axis_magnitude;

    Tensor result;
    result.xx() = (q0 * q0) + (q1 * q1) - (q2 * q2) - (q3 * q3);
    result.xy() = 2.0_rt * (q1 * q2 - q0 * q3);
    result.xz() = 2.0_rt * (q0 * q2 + q1 * q3);

    result.yx() = 2.0_rt * (q1 * q2 + q0 * q3);
    result.yy() = (q0 * q0) - (q1 * q1) + (q2 * q2) - (q3 * q3);
    result.yz() = 2.0_rt * (q2 * q3 - q0 * q1);

    result.zx() = 2.0_rt * (q1 * q3 - q0 * q2);
    result.zy() = 2.0_rt * (q0 * q1 + q2 * q3);
    result.zz() = (q0 * q0) - (q1 * q1) - (q2 * q2) + (q3 * q3);

    return result;
}

/** Construct Rz(yaw) Ry(pitch) Rx(roll), with angles in degrees. */
inline Quaternion from_roll_pitch_yaw(const Vector& angles)
{
    // The existing vs rotation tensors use the opposite sign from the
    // conventional active quaternion formula.
    const amrex::Real roll = -0.5_rt * utils::radians(angles.x());
    const amrex::Real pitch = -0.5_rt * utils::radians(angles.y());
    const amrex::Real yaw = -0.5_rt * utils::radians(angles.z());
    const amrex::Real cr = std::cos(roll);
    const amrex::Real sr = std::sin(roll);
    const amrex::Real cp = std::cos(pitch);
    const amrex::Real sp = std::sin(pitch);
    const amrex::Real cy = std::cos(yaw);
    const amrex::Real sy = std::sin(yaw);
    return normalized(
        {cr * cp * cy + sr * sp * sy, sr * cp * cy - cr * sp * sy,
         cr * sp * cy + sr * cp * sy, cr * cp * sy - sr * sp * cy});
}

inline Quaternion slerp(Quaternion a, Quaternion b, const amrex::Real fraction)
{
    a = normalized(a);
    b = normalized(b);
    amrex::Real cosine = dot(a, b);
    if (cosine < 0.0_rt) {
        b = {-b.w, -b.x, -b.y, -b.z};
        cosine = -cosine;
    }
    // Scale the near-parallel cutoff with the configured Real precision.
    if (cosine > 1.0_rt - std::sqrt(constants::EPS)) {
        return normalized(
            {a.w + fraction * (b.w - a.w), a.x + fraction * (b.x - a.x),
             a.y + fraction * (b.y - a.y), a.z + fraction * (b.z - a.z)});
    }
    const amrex::Real angle = std::acos(std::clamp(cosine, -1.0_rt, 1.0_rt));
    const amrex::Real denom = std::sin(angle);
    const amrex::Real wa = std::sin((1.0_rt - fraction) * angle) / denom;
    const amrex::Real wb = std::sin(fraction * angle) / denom;
    return {
        wa * a.w + wb * b.w, wa * a.x + wb * b.x, wa * a.y + wb * b.y,
        wa * a.z + wb * b.z};
}

} // namespace kynema_sgf::vs

#endif
