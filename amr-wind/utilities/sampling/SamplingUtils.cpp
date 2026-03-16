#include <algorithm>

#include "SamplingUtils.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::sampling::sampling_utils {

vs::Vector reflect(vs::Vector line, vs::Vector vec)
{
    vs::Tensor ref(
        1 - (2 * line.x() * line.x()), -2 * line.x() * line.y(),
        -2 * line.x() * line.z(), -2 * line.y() * line.x(),
        1 - (2 * line.y() * line.y()), -2 * line.y() * line.z(),
        -2 * line.z() * line.x(), -2 * line.z() * line.y(),
        1 - (2 * line.z() * line.z()));

    return vec & ref;
}

vs::Vector rotate_euler_vec(vs::Vector axis, amrex::Real angle, vs::Vector vec)
{
    axis.normalize();
    const auto RotMat = vs::quaternion(axis, angle);
    return vec & RotMat;
}

vs::Vector
rotate_euler_vector(vs::Vector& axis, amrex::Real& angle, vs::Vector& vec)
{
    axis.normalize();
    const auto RotMat = vs::quaternion(axis, angle);
    return vec & RotMat;
}

vs::Vector rotation(const vs::Vector& angles, const vs::Vector& data)
{
    const vs::Tensor rotMatrix =
        vs::xrot(angles.x()) & vs::yrot(angles.y()) & vs::zrot(angles.z());
    return data & rotMatrix;
}

vs::Vector canon_rotator(const vs::Vector& angles, const vs::Vector& data)
{
    // Note: Angles in degrees as per "vs" class
    const vs::Tensor rotMatrix =
        vs::xrot(angles.x()) & vs::yrot(angles.y()) & vs::zrot(angles.z());
    return data & rotMatrix;
}

vs::Tensor unit_projection_matrix(const vs::Vector& a)
{
    return {a[0] * a[0], a[0] * a[1], a[0] * a[2], a[0] * a[1], a[1] * a[1],
            a[1] * a[2], a[0] * a[2], a[1] * a[2], a[2] * a[2]};
}

vs::Tensor rotation_matrix(vs::Vector dst, vs::Vector src)
{
    auto vmat = skew_cross(dst, src);
    const auto ang = dst & src;

    const amrex::Real small =
        std::numeric_limits<amrex::Real>::epsilon() * 1.0e2_rt * vs::mag(dst);
    if (std::abs(1 + ang) < small) {
        return scale(vs::Tensor::identity(), -1);
    }
    return vs::Tensor::identity() + vmat +
           scale((vmat & vmat), 1.0_rt / (1 + ang));
}

vs::Tensor skew_cross(vs::Vector a, vs::Vector b)
{
    auto cross = b ^ a;
    return {0,         -cross[2], cross[1], cross[2], 0,
            -cross[0], -cross[1], cross[0], 0};
}

vs::Tensor scale(vs::Tensor v, amrex::Real a)
{
    vs::Tensor vnew;
    for (int j = 0; j < 9; ++j) {
        vnew[j] = a * v[j];
    }
    return vnew;
}

void spherical_cap_quadrature(
    amrex::Real gammav,
    int ntheta,
    std::vector<amrex::Real> abscissae1D,
    std::vector<amrex::Real> weights1D,
    std::vector<vs::Vector>& rays,
    std::vector<amrex::Real>& weights)
{
    // tensor-product quadrature. theta is uniform on [0, 2 pi) since it's a
    // circle and phi is a standard 1D quadrature of choice, mapped to
    // [cos(gamma), 1] with a variable transformation, following along from
    // DOI 10.1007/s10444-011-9187-2
    std::ranges::transform(
        abscissae1D, abscissae1D.begin(), [gammav](amrex::Real s) {
            return (0.5_rt * (1.0_rt - std::cos(gammav)) * (-s)) +
                   (0.5_rt * (1.0_rt + std::cos(gammav)));
        });
    std::ranges::transform(
        weights1D, weights1D.begin(), [gammav](amrex::Real w) {
            return w * 0.5_rt * (1.0_rt - std::cos(gammav));
        });

    const int nphi = int(abscissae1D.size());

    // avoid ntheta multiplicity at the center
    rays[0] = vs::Vector(0, 0, 1);
    weights[0] = (2.0_rt * std::numbers::pi_v<amrex::Real> * weights1D[0]);

    const auto theta_weight = 2.0_rt * std::numbers::pi_v<amrex::Real> / ntheta;
    for (int j = 1; j < nphi; ++j) {
        const auto tau = abscissae1D[j];
        for (int i = 0; i < ntheta; ++i) {
            int r_idx = i + ((j - 1) * ntheta) + 1;
            const amrex::Real theta =
                (2.0_rt * std::numbers::pi_v<amrex::Real> / ntheta) * i;
            const auto xr = std::sqrt(1.0_rt - (tau * tau)) * std::cos(theta);
            const auto yr = std::sqrt(1.0_rt - (tau * tau)) * std::sin(theta);
            const auto zr = tau;
            auto ray = vs::Vector(xr, yr, zr);
            ray.normalize();
            rays[r_idx] = ray;
            const auto phi_weight = weights1D[j];
            weights[r_idx] = theta_weight * phi_weight;
        }
    }
}

void spherical_cap_truncated_normal(
    amrex::Real gammav,
    int ntheta,
    NormalRule rule,
    std::vector<vs::Vector>& rays,
    std::vector<amrex::Real>& weights)
{
    auto [xlocs, xweights] = truncated_normal_rule(rule);
    // want the center of the truncated normal distribution at the pole of the
    // cap -> -1 . Weights are already for a [-1,1] range from the generator
    std::ranges::transform(
        xlocs, xlocs.begin(), [](amrex::Real x) { return (2 * x) - 1; });
    // half range to start, then mapped back to [-1,1]
    std::ranges::transform(
        xweights, xweights.begin(), [](amrex::Real w) { return 4 * w; });

    spherical_cap_quadrature(gammav, ntheta, xlocs, xweights, rays, weights);
}

std::pair<std::vector<amrex::Real>, std::vector<amrex::Real>>
truncated_normal_rule(NormalRule rule)
{
    // from the "truncated normal quadrature" .python code
    switch (rule) {
    case NormalRule::SIGMA1:
        return {
            {0, 0.3436121352489559_rt, 0.6473220334297102_rt,
             0.8706217027202311_rt, 0.9816974860670653_rt},
            {0.2046394226558322_rt / 2.0_rt, 0.1820146209511494_rt,
             0.128027596313765_rt, 0.06821017522834351_rt,
             0.01942789617882675_rt}};
    case NormalRule::SIGMA2:
        return {
            {0, 0.2959590846054073_rt, 0.5735693238435292_rt,
             0.8074757570903542_rt, 0.9607561326630086_rt},
            {0.249758577881117_rt / 2.0_rt, 0.2035976917174024_rt,
             0.1129523637830892_rt, 0.04552496709664563_rt,
             0.013045688462303995_rt}};
    case NormalRule::SIGMA3:
        return {
            {0, 0.2661790968493327_rt, 0.5263305051027921_rt,
             0.7664900509748058_rt, 0.9477581057921652_rt},
            {0.3203929665957703_rt / 2.0_rt, 0.2307493381206665_rt,
             0.08754316928625956_rt, 0.01882073900490792_rt,
             0.002690270290280566_rt}};
    case NormalRule::HALFPOWER:
        return {
            {0, 0.315493297131259_rt, 0.6016636608468_rt, 0.8282105821126121_rt,
             0.9662550592631028_rt},
            {0.197723576944154_rt / 2.0_rt, 0.1761766447490471_rt,
             0.1255723775152601_rt, 0.07163437433902098_rt,
             0.02775481492459504_rt}};
    default: {
        throw std::runtime_error(
            "Only implemented 1-3, halfpower for truncated normal");
        return {};
    }
    }
}

} // namespace amr_wind::sampling::sampling_utils
