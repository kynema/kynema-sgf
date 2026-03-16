#ifndef SAMPLINGUTILS_H
#define SAMPLINGUTILS_H

#include <string>
#include <cstdint>
#include "amr-wind/core/vs/vector_space.H"
#include "amr-wind/utilities/tensor_ops.H"
#include "AMReX_Vector.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::sampling::sampling_utils {

/** Utility functions available to advanced sampling techniques
 *
 *  \param key Name of the PDE
 *  \param terms List of identifiers for the source terms to be added
 */
// void inject_turbulence_src_terms(
//     const std::string& key, const amrex::Vector<std::string>& terms);

enum class NormalRule : std::uint8_t { SIGMA1, SIGMA2, SIGMA3, HALFPOWER };

vs::Vector reflect(vs::Vector line, vs::Vector vec);
vs::Vector rotate_euler_vec(vs::Vector axis, amrex::Real angle, vs::Vector vec);
vs::Vector
rotate_euler_vector(vs::Vector& axis, amrex::Real& angle, vs::Vector& vec);
vs::Vector rotation(const vs::Vector& angles, const vs::Vector& data);
vs::Vector canon_rotator(const vs::Vector& angles, const vs::Vector& data);
vs::Tensor unit_projection_matrix(const vs::Vector& a);
vs::Tensor rotation_matrix(vs::Vector dst, vs::Vector src);
vs::Tensor skew_cross(vs::Vector a, vs::Vector b);
vs::Tensor scale(vs::Tensor v, amrex::Real a);

void spherical_cap_quadrature(
    amrex::Real gammav,
    int ntheta,
    std::vector<amrex::Real> abscissae1D,
    std::vector<amrex::Real> weights1D,
    std::vector<vs::Vector>& rays,
    std::vector<amrex::Real>& weights);

void spherical_cap_truncated_normal(
    amrex::Real gammav,
    int ntheta,
    NormalRule rule,
    std::vector<vs::Vector>& rays,
    std::vector<amrex::Real>& weights);

std::pair<std::vector<amrex::Real>, std::vector<amrex::Real>>
truncated_normal_rule(NormalRule rule);

} // namespace amr_wind::sampling::sampling_utils

#endif /* SAMPLINGUTILS_H */
