#ifndef FLLCOP_H
#define FLLCOP_H

#include <numbers>
#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/FLLC.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator {

/** \brief This struct will operate on a blade/wing.
 *  The velocity from the simulation is corrected using the
 *  Filtered Lifting Line Theory Correction from Martinez and Meneveau 2019
 *
 */

struct FLLCOp
{
    void operator()(ComponentView& data, FLLCData& fllc)
    {
        if (!fllc.initialized) {
            return;
        }
        switch (fllc.correction_type) {
        case FLLCType::ConstantChord: {
            constant_chord(data, fllc);
            break;
        }
        case FLLCType::VariableChord: {
            variable_chord(data, fllc);
            break;
        }
        default:
            // this should be unreachable
            throw std::runtime_error("Invalid FLLC type");
        }
    }

    static void constant_chord(ComponentView& data, FLLCData& fllc)
    {
        namespace interp = ::amr_wind::interp;

        const int npts = static_cast<int>(fllc.correction_velocity.size());
        auto& du = fllc.correction_velocity;

        auto& u_les = fllc.les_velocity;
        auto& u_opt = fllc.optimal_velocity;
        auto& G = fllc.lift;
        auto& dG = fllc.grad_lift;
        VecSlice u_force_pnt_slice =
            ::amr_wind::utils::slice(fllc.force_point_velocity, 0);
        vs::Vector* vel_ptr = data.vel.data();

        /**
         * Step 0
         * Linear interpolate velocity if force and vel points are mismatched
         */
        if (fllc.different_sizes) {
            interp::linear(
                fllc.span_distance_vel, data.vel, fllc.span_distance_force,
                u_force_pnt_slice);
            vel_ptr = u_force_pnt_slice.data();
        }

        /**
         * Step 1
         * Compute the lift force distribution (G)
         * Compute equation 5.3_rt from Martinez-Tossas and Meneveau 2019
         */
        for (int ip = 0; ip < npts; ++ip) {
            const auto force = data.force[ip];
            const auto vel = data.vel_rel[ip];
            const auto dr = fllc.dr[ip];
            const auto vmag = amrex::max<amrex::Real>(
                vs::mag(vel), vs::DTraits<amrex::Real>::eps());
            const auto vmag2 = vmag * vmag;

            const auto fv = force & vel;

            G[ip] = (force - vel * fv / vmag2) / dr;
        }

        /**
         * Step 2
         * Compute gradient of the lift force distribution (\Delta G)
         * Compute equations 5.4_rt and 5.5_rt from Martinez-Tossas and Meneveau
         * 2019
         */
        dG[0] = G[0];
        dG[npts - 1] = -1 * G[npts - 1];
        for (int ip = 1; ip < npts - 1; ++ip) {
            dG[ip] = 0.5_rt * (G[ip + 1] - G[ip - 1]);
        }

        /**
         * Step 3
         * Compute the induced velocities
         * Compute equations 5.6_rt and 5.7_rt from Martinez-Tossas and Meneveau
         * 2019
         */
        for (int ip = 0; ip < npts; ++ip) {

            const auto eps_les = fllc.epsilon;
            const auto eps_opt = fllc.optimal_epsilon[ip];

            for (int jp = 0; jp < npts; ++jp) {

                if (ip == jp) {
                    continue;
                }

                const auto r_dis = vs::mag(data.vel_pos[ip] - data.vel_pos[jp]);
                const auto& vel = data.vel_rel[jp];
                const auto vmag = amrex::max<amrex::Real>(
                    vs::mag(vel), vs::DTraits<amrex::Real>::eps());

                auto coefficient =
                    1.0_rt /
                    (-4.0_rt * std::numbers::pi_v<amrex::Real> * r_dis * vmag);
                const auto cLes =
                    1.0_rt - std::exp(-r_dis * r_dis / (eps_les * eps_les));
                const auto cOpt =
                    1.0_rt - std::exp(-r_dis * r_dis / (eps_opt * eps_opt));

                // The sign of the induced velocity depends on which side of the
                // blade we are on
                if (ip < jp) {
                    coefficient *= -1.0_rt;
                }

                u_les[ip] = u_les[ip] - dG[jp] * coefficient * cLes;
                u_opt[ip] = u_opt[ip] - dG[jp] * coefficient * cOpt;
            }

            // Relaxation to compute the induced velocity
            const auto f = fllc.relaxation_factor;
            du[ip] = (1.0_rt - f) * du[ip] + f * (u_opt[ip] - u_les[ip]);
        }

        /**
         * Step 4
         * Compute the induced velocity difference
         * Compute equation 5.8_rt from Martinez-Tossas and Meneveau 2019
         */
        for (int ip = 0; ip < npts; ++ip) {
            vel_ptr[ip] = vel_ptr[ip] + du[ip];
            u_les[ip] *= 0.0_rt;
            u_opt[ip] *= 0.0_rt;
        }

        /**
         * Step 5
         * Interpolate back to vel positions if necessary
         */
        if (fllc.different_sizes) {
            interp::linear(
                fllc.span_distance_force, u_force_pnt_slice,
                fllc.span_distance_vel, data.vel);
        }
    }

    static void variable_chord(ComponentView& data, FLLCData& fllc)
    {
        namespace interp = ::amr_wind::interp;

        const int npts = static_cast<int>(fllc.correction_velocity.size());
        auto& du = fllc.correction_velocity;

        auto& nonuniform_r = fllc.nonuniform_r;
        auto& u_les = fllc.les_velocity;
        auto& u_opt = fllc.optimal_velocity;
        auto& G = fllc.lift;
        auto& nonuniform_G = fllc.nonuniform_lift;
        VecSlice u_force_pnt_slice =
            ::amr_wind::utils::slice(fllc.force_point_velocity, 0);
        vs::Vector* vel_ptr = data.vel.data();

        /**
         * Step 0
         * Linear interpolate velocity if force and vel points are mismatched
         */
        if (fllc.different_sizes) {
            interp::linear(
                fllc.span_distance_vel, data.vel, fllc.span_distance_force,
                u_force_pnt_slice);
            vel_ptr = u_force_pnt_slice.data();
        }

        /**
         * Step 1
         * Compute the lift force distribution (G)
         */
        for (int ip = 0; ip < npts; ++ip) {
            const auto force = data.force[ip];
            const auto vel = data.vel_rel[ip];
            const auto dr = fllc.dr[ip];
            const auto vmag = amrex::max<amrex::Real>(
                vs::mag(vel), vs::DTraits<amrex::Real>::eps());
            const auto vmag2 = vmag * vmag;
            fllc.vel_rel[ip] = vel;
            const auto fv = force & vel;

            G[ip] = (force - vel * fv / vmag2) / dr;
        }

        // Get the non-uniform distribution of force
        if (fllc.nonuniform) {
            interp::linear(fllc.r, G, nonuniform_r, nonuniform_G);
            interp::linear(
                fllc.r, fllc.optimal_epsilon, nonuniform_r,
                fllc.nonuniform_optimal_epsilon);
            interp::linear(
                fllc.r, fllc.vel_rel, nonuniform_r, fllc.nonuniform_vel_rel);
        }

        const auto optimal_epsilon = (fllc.nonuniform)
                                         ? fllc.nonuniform_optimal_epsilon
                                         : fllc.optimal_epsilon;
        const auto dr = (fllc.nonuniform) ? fllc.nonuniform_dr : fllc.dr;
        const auto vel_rel =
            (fllc.nonuniform) ? fllc.nonuniform_vel_rel : fllc.vel_rel;
        const auto G_new = (fllc.nonuniform) ? nonuniform_G : G;

        /**
         * Step 2
         * Compute the induced velocities
         */
        for (int ip = 0; ip < npts; ++ip) {

            const auto coefficient =
                1.0_rt / (2.0_rt * std::numbers::pi_v<amrex::Real>);

            for (int jp = 0; jp < dr.size(); ++jp) {

                const auto eps_les = fllc.epsilon;
                const auto eps_opt = optimal_epsilon[jp];
                const auto eps_les2 = eps_les * eps_les;
                const auto eps_opt2 = eps_opt * eps_opt;
                const auto& vel = vel_rel[jp];
                const auto vmag = amrex::max<amrex::Real>(
                    vs::mag(vel), vs::DTraits<amrex::Real>::eps());
                auto k_les = 0.0_rt;
                auto k_opt = 0.0_rt;

                const auto r_dis =
                    (fllc.nonuniform)
                        ? std::abs(fllc.r[ip] - fllc.nonuniform_r[jp])
                        : vs::mag(data.vel_pos[ip] - data.vel_pos[jp]);

                if (r_dis == 0.) {

                    k_les = 0.5_rt / (eps_les2);
                    k_opt = 0.5_rt / (eps_opt2);

                }

                else {

                    auto exp_les = std::exp(-r_dis * r_dis / eps_les2);
                    auto exp_opt = std::exp(-r_dis * r_dis / eps_opt2);

                    k_les = (exp_les / eps_les2) +
                            (1.0_rt / (2.0_rt * r_dis * r_dis) *
                             (exp_les - 1.0_rt));
                    k_opt = (exp_opt / eps_opt2) +
                            (1.0_rt / (2.0_rt * r_dis * r_dis) *
                             (exp_opt - 1.0_rt));
                }

                u_les[ip] =
                    u_les[ip] + coefficient * G_new[jp] / vmag * k_les * dr[jp];
                u_opt[ip] =
                    u_opt[ip] + coefficient * G_new[jp] / vmag * k_opt * dr[jp];
            }

            // Relaxation to compute the induced velocity
            const auto f = fllc.relaxation_factor;
            du[ip] = (1.0_rt - f) * du[ip] + f * (u_opt[ip] - u_les[ip]);
        }

        /**
         * Step 3
         * Compute the induced velocity difference
         */
        for (int ip = 0; ip < npts; ++ip) {
            vel_ptr[ip] = vel_ptr[ip] + du[ip];
            u_les[ip] *= 0.0_rt;
            u_opt[ip] *= 0.0_rt;
        }

        /**
         * Step 4
         * Interpolate back to vel positions if necessary
         */
        if (fllc.different_sizes) {
            interp::linear(
                fllc.span_distance_force, u_force_pnt_slice,
                fllc.span_distance_vel, data.vel);
        }
    }
};

} // namespace amr_wind::actuator

#endif /* FLLCOP_H */
