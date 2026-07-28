#ifndef VS_QUATERNION_H
#define VS_QUATERNION_H

#include "src/core/vs/tensor.H"
#include "src/utilities/trig_ops.H"

#include <cmath>

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

    //! Normalize this quaternion in place
    Quaternion& normalize();

    //! Return a normalized copy without modifying this quaternion
    [[nodiscard]] Quaternion normalized() const;

    //! Return the quaternion conjugate
    [[nodiscard]] Quaternion conjugate() const;
};

amrex::Real dot(const Quaternion& a, const Quaternion& b);

//! Rotation composition: apply rhs first, followed by lhs
Quaternion operator*(const Quaternion& lhs, const Quaternion& rhs);

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

Tensor tensor(const Quaternion& input);

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
Quaternion from_roll_pitch_yaw(const Vector& angles);

Quaternion slerp(Quaternion a, Quaternion b, amrex::Real fraction);

} // namespace kynema_sgf::vs

#endif
