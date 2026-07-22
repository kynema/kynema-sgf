#ifndef TERRAINDRAGDAMPING_H
#define TERRAINDRAGDAMPING_H

#include <cmath>
#include <numbers>

#include "AMReX_Geometry.H"
#include "AMReX_ParmParse.H"
#include "AMReX_REAL.H"
#include "src/core/Field.H"

namespace kynema_sgf::terraindrag::damping_legacy {

inline bool has_damping_inputs(const amrex::ParmParse& pp)
{
    return pp.contains("damp_east_slope") || pp.contains("damp_east_full") ||
           pp.contains("damp_west_slope") || pp.contains("damp_west_full") ||
           pp.contains("damp_north_slope") ||
           pp.contains("damp_north_full") ||
           pp.contains("damp_south_slope") ||
           pp.contains("damp_south_full") ||
           pp.contains("horizontal_time_scale") ||
           pp.contains("horizontal_abl_height") ||
           pp.contains("horizontal_slope_end") ||
           pp.contains("vertical_slope") || pp.contains("vertical_full");
}

inline void query_damping_inputs(
    amrex::ParmParse& pp,
    amrex::Real& damp_east_slope,
    amrex::Real& damp_east_full,
    amrex::Real& damp_west_slope,
    amrex::Real& damp_west_full,
    amrex::Real& damp_north_slope,
    amrex::Real& damp_north_full,
    amrex::Real& damp_south_slope,
    amrex::Real& damp_south_full,
    amrex::Real& horizontal_tau,
    amrex::Real& horizontal_abl_height,
    amrex::Real& horizontal_slope_end,
    amrex::Real& vertical_slope,
    amrex::Real& vertical_full)
{
    pp.query("damp_east_slope", damp_east_slope);
    pp.query("damp_east_full", damp_east_full);
    pp.query("damp_west_slope", damp_west_slope);
    pp.query("damp_west_full", damp_west_full);
    pp.query("damp_north_slope", damp_north_slope);
    pp.query("damp_north_full", damp_north_full);
    pp.query("damp_south_slope", damp_south_slope);
    pp.query("damp_south_full", damp_south_full);
    pp.query("horizontal_time_scale", horizontal_tau);
    pp.query("horizontal_abl_height", horizontal_abl_height);
    pp.query("horizontal_slope_end", horizontal_slope_end);
    pp.query("vertical_slope", vertical_slope);
    pp.query("vertical_full", vertical_full);
}

inline void compute_terrain_damping(
    Field& terrain_damping,
    int level,
    const amrex::Geometry& geom,
    amrex::Real damp_east_slope,
    amrex::Real damp_east_full,
    amrex::Real damp_west_slope,
    amrex::Real damp_west_full,
    amrex::Real damp_north_slope,
    amrex::Real damp_north_full,
    amrex::Real damp_south_slope,
    amrex::Real damp_south_full,
    amrex::Real horizontal_tau,
    amrex::Real horizontal_abl_height,
    amrex::Real horizontal_slope_end,
    amrex::Real vertical_slope,
    amrex::Real vertical_full)
{
    using namespace amrex::literals;

    auto& damping = terrain_damping(level);
    auto level_damping = damping.arrays();

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();

    const amrex::Real damping_east_start =
        prob_hi[0] - (damp_east_full + damp_east_slope);
    const amrex::Real damping_east_end = prob_hi[0] - damp_east_full;
    const amrex::Real damping_west_start =
        prob_lo[0] + (damp_west_full + damp_west_slope);
    const amrex::Real damping_west_end = prob_lo[0] + damp_west_full;
    const amrex::Real damping_north_start =
        prob_hi[1] - (damp_north_full + damp_north_slope);
    const amrex::Real damping_north_end = prob_hi[1] - damp_north_full;
    const amrex::Real damping_south_start =
        prob_lo[1] + (damp_south_full + damp_south_slope);
    const amrex::Real damping_south_end = prob_lo[1] + damp_south_full;

    amrex::ParallelFor(
        damping, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            amrex::Real horizontal_coeff_east = 0.0_rt;
            amrex::Real horizontal_coeff_north = 0.0_rt;
            amrex::Real horizontal_coeff_west = 0.0_rt;
            amrex::Real horizontal_coeff_south = 0.0_rt;
            amrex::Real vertical_coeff = 0.0_rt;
            const amrex::Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
            const amrex::Real y = prob_lo[1] + (j + 0.5_rt) * dx[1];
            const amrex::Real z = prob_lo[2] + (k + 0.5_rt) * dx[2];

            if (x < damping_east_start) {
                horizontal_coeff_east = 0.0_rt;
            } else if (x >= damping_east_end) {
                horizontal_coeff_east = 1.0_rt;
            } else {
                const amrex::Real term = std::sin(
                    std::numbers::pi_v<amrex::Real> * 0.5_rt *
                    (x - damping_east_start) /
                    (damping_east_end - damping_east_start));
                horizontal_coeff_east = term * term;
            }

            if (x > damping_west_start) {
                horizontal_coeff_west = 0.0_rt;
            } else if (x <= damping_west_end) {
                horizontal_coeff_west = 1.0_rt;
            } else {
                const amrex::Real term = std::sin(
                    std::numbers::pi_v<amrex::Real> * 0.5_rt *
                    (x - damping_west_start) /
                    (damping_west_end - damping_west_start));
                horizontal_coeff_west = term * term;
            }

            if (y < damping_north_start) {
                horizontal_coeff_north = 0.0_rt;
            } else if (y >= damping_north_end) {
                horizontal_coeff_north = 1.0_rt;
            } else {
                const amrex::Real term = std::sin(
                    std::numbers::pi_v<amrex::Real> * 0.5_rt *
                    (y - damping_north_start) /
                    (damping_north_end - damping_north_start));
                horizontal_coeff_north = term * term;
            }

            if (y > damping_south_start) {
                horizontal_coeff_south = 0.0_rt;
            } else if (y <= damping_south_end) {
                horizontal_coeff_south = 1.0_rt;
            } else {
                const amrex::Real term = std::sin(
                    std::numbers::pi_v<amrex::Real> * 0.5_rt *
                    (y - damping_south_start) /
                    (damping_south_end - damping_south_start));
                horizontal_coeff_south = term * term;
            }

            if (z <= horizontal_abl_height) {
                vertical_coeff = 0.0_rt;
            } else if (z > horizontal_slope_end) {
                vertical_coeff = 1.0_rt;
            } else {
                const amrex::Real term = std::sin(
                    std::numbers::pi_v<amrex::Real> * 0.5_rt *
                    (z - horizontal_abl_height) /
                    (horizontal_slope_end - horizontal_abl_height));
                vertical_coeff = term * term;
            }

            level_damping[nbx](i, j, k, 0) =
                vertical_coeff *
                (horizontal_coeff_east + horizontal_coeff_north +
                 horizontal_coeff_west + horizontal_coeff_south);

            if (z <= vertical_slope) {
                vertical_coeff = 0.0_rt;
            } else if (z > vertical_full) {
                vertical_coeff = 1.0_rt;
            } else {
                const amrex::Real term = std::sin(
                    std::numbers::pi_v<amrex::Real> * 0.5_rt *
                    (z - vertical_slope) /
                    (vertical_full - vertical_slope + 1e-15_rt));
                vertical_coeff = term * term;
            }

            level_damping[nbx](i, j, k, 0) =
                amrex::min<amrex::Real>(
                    vertical_coeff + level_damping[nbx](i, j, k, 0), 1.0_rt) /
                horizontal_tau;
        });
}

} // namespace kynema_sgf::terraindrag::damping_legacy

#endif
