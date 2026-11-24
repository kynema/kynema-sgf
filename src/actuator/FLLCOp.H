#ifndef FLLCOP_H
#define FLLCOP_H

#include <cmath>
#include <numbers>
#include "src/actuator/actuator_types.H"
#include "src/actuator/FLLC.H"
#include "src/utilities/linear_interpolation.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::actuator {

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
        const int npts = static_cast<int>(fllc.correction_velocity.size());
        auto& du = fllc.correction_velocity;

        auto& u_les = fllc.les_velocity;
        auto& u_opt = fllc.optimal_velocity;
        auto& G = fllc.lift;
        auto& dG = fllc.grad_lift;

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

            // The force is spread with the largest of epsilon, epsilon_min
            // (turbines; fllc.epsilon is the larger of the two) and the
            // optimal epsilon, so that is the width the simulation resolves
            const auto eps_opt = fllc.optimal_epsilon[ip];
            const auto eps_les = amrex::max<amrex::Real>(fllc.epsilon, eps_opt);

            for (int jp = 0; jp < npts; ++jp) {

                if (ip == jp) {
                    continue;
                }

                const auto r_dis = vs::mag(data.pos[ip] - data.pos[jp]);
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
        add_correction(data, fllc);
        for (int ip = 0; ip < npts; ++ip) {
            u_les[ip] *= 0.0_rt;
            u_opt[ip] *= 0.0_rt;
        }
    }

    /** Add the correction to the sampled velocity
     *
     *  The correction is known at the force points and the velocity is
     *  sampled at the velocity points. Where the two sets differ, the
     *  correction is interpolated along the span to the velocity points; the
     *  sampled velocity itself is never interpolated.
     */
    static void add_correction(ComponentView& data, const FLLCData& fllc)
    {
        namespace interp = ::kynema_sgf::interp;

        const auto& du = fllc.correction_velocity;
        const int npts_vel = static_cast<int>(data.vel.size());
        for (int ip = 0; ip < npts_vel; ++ip) {
            if (fllc.mismatched_points) {
                const auto du_vel = interp::linear(
                    fllc.span_distance_force, du, fllc.span_distance_vel[ip]);
                data.vel[ip] = data.vel[ip] + du_vel;
            } else {
                data.vel[ip] = data.vel[ip] + du[ip];
            }
        }
    }

    /** Kernel of the velocity induced by a Gaussian-filtered lifting line
     *
     *  exp(-r^2 / eps^2) / eps^2 + (exp(-r^2 / eps^2) - 1) / (2 r^2), with
     *  the limit 1 / (2 eps^2) at r = 0. The second term is evaluated with
     *  expm1: written as exp(x) - 1 it loses all its digits when r is much
     *  smaller than epsilon, and the kernel is then twice its limit.
     *
     *  \param r_dis2 Square of the distance between the two points
     *  \param eps2 Square of epsilon
     */
    static amrex::Real kernel(const amrex::Real r_dis2, const amrex::Real eps2)
    {
        const auto x = -r_dis2 / eps2;
        if (x == 0.0_rt) {
            return 0.5_rt / eps2;
        }
        return (std::exp(x) / eps2) + (std::expm1(x) / (2.0_rt * r_dis2));
    }

    static void variable_chord(ComponentView& data, FLLCData& fllc)
    {
        namespace interp = ::kynema_sgf::interp;

        const int npts = static_cast<int>(fllc.correction_velocity.size());
        auto& du = fllc.correction_velocity;

        auto& nonuniform_r = fllc.nonuniform_r;
        auto& u_les = fllc.les_velocity;
        auto& u_opt = fllc.optimal_velocity;
        auto& G = fllc.lift;
        auto& nonuniform_G = fllc.nonuniform_lift;

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

        const auto& optimal_epsilon = (fllc.nonuniform)
                                          ? fllc.nonuniform_optimal_epsilon
                                          : fllc.optimal_epsilon;
        const auto& dr = (fllc.nonuniform) ? fllc.nonuniform_dr : fllc.dr;
        const auto& vel_rel =
            (fllc.nonuniform) ? fllc.nonuniform_vel_rel : fllc.vel_rel;
        const auto& G_new = (fllc.nonuniform) ? nonuniform_G : G;

        /**
         * Step 2
         * Compute the induced velocities
         */
        for (int ip = 0; ip < npts; ++ip) {

            const auto coefficient =
                1.0_rt / (2.0_rt * std::numbers::pi_v<amrex::Real>);

            for (int jp = 0; jp < dr.size(); ++jp) {

                // The force is spread with the largest of epsilon,
                // epsilon_min (turbines; fllc.epsilon is the larger of the
                // two) and the optimal epsilon, so that is the width the
                // simulation resolves
                const auto eps_opt = optimal_epsilon[jp];
                const auto eps_les =
                    amrex::max<amrex::Real>(fllc.epsilon, eps_opt);
                const auto eps_les2 = eps_les * eps_les;
                const auto eps_opt2 = eps_opt * eps_opt;
                const auto& vel = vel_rel[jp];
                const auto vmag = amrex::max<amrex::Real>(
                    vs::mag(vel), vs::DTraits<amrex::Real>::eps());

                const auto r_dis =
                    (fllc.nonuniform)
                        ? std::abs(fllc.r[ip] - fllc.nonuniform_r[jp])
                        : vs::mag(data.pos[ip] - data.pos[jp]);
                const auto r_dis2 = r_dis * r_dis;
                const auto k_les = kernel(r_dis2, eps_les2);
                const auto k_opt = kernel(r_dis2, eps_opt2);

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
        add_correction(data, fllc);
        for (int ip = 0; ip < npts; ++ip) {
            u_les[ip] *= 0.0_rt;
            u_opt[ip] *= 0.0_rt;
        }
    }
};

} // namespace kynema_sgf::actuator

#endif /* FLLCOP_H */
