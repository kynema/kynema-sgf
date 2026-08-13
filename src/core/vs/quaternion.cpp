#include "src/core/vs/quaternion.H"
#include "src/utilities/constants.H"

#include <algorithm>
#include <cmath>

namespace kynema_sgf::vs {

using namespace amrex::literals;

Quaternion& Quaternion::normalize()
{
    const amrex::Real norm = std::sqrt(w * w + x * x + y * y + z * z);
    if (norm <= constants::EPS) {
        amrex::Abort("Cannot normalize a zero quaternion");
    }
    w /= norm;
    x /= norm;
    y /= norm;
    z /= norm;
    return *this;
}

Quaternion Quaternion::normalized() const
{
    Quaternion result = *this;
    result.normalize();
    return result;
}

Quaternion Quaternion::conjugate() const { return {w, -x, -y, -z}; }

amrex::Real dot(const Quaternion& a, const Quaternion& b)
{
    return a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
}

Quaternion operator*(const Quaternion& lhs, const Quaternion& rhs)
{
    return {
        lhs.w * rhs.w - lhs.x * rhs.x - lhs.y * rhs.y - lhs.z * rhs.z,
        lhs.w * rhs.x + lhs.x * rhs.w + lhs.y * rhs.z - lhs.z * rhs.y,
        lhs.w * rhs.y - lhs.x * rhs.z + lhs.y * rhs.w + lhs.z * rhs.x,
        lhs.w * rhs.z + lhs.x * rhs.y - lhs.y * rhs.x + lhs.z * rhs.w};
}

Tensor tensor(const Quaternion& input)
{
    return tensor_unit(input.normalized());
}

Quaternion from_roll_pitch_yaw(const Vector& angles)
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
    return Quaternion{
        cr * cp * cy + sr * sp * sy, sr * cp * cy - cr * sp * sy,
        cr * sp * cy + sr * cp * sy, cr * cp * sy - sr * sp * cy}
        .normalized();
}

Quaternion slerp(Quaternion a, Quaternion b, const amrex::Real fraction)
{
    a.normalize();
    b.normalize();
    amrex::Real cosine = dot(a, b);
    // Flip one representation when needed so interpolation follows the
    // shortest path between orientations.
    if (cosine < 0.0_rt) {
        b = {-b.w, -b.x, -b.y, -b.z};
        cosine = -cosine;
    }
    // Scale the near-parallel cutoff with the configured Real precision.
    if (cosine > 1.0_rt - std::sqrt(constants::EPS)) {
        return Quaternion{
            a.w + fraction * (b.w - a.w), a.x + fraction * (b.x - a.x),
            a.y + fraction * (b.y - a.y), a.z + fraction * (b.z - a.z)}
            .normalized();
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
