#include <cmath>
#include <limits>
#include <sstream>
#include "src/actuator/FLLC.H"
#include "src/utilities/linear_interpolation.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::actuator {

void fllc_init(
    FLLCData& data,
    const ComponentView& view,
    const amrex::Real eps_chord,
    const std::string& label,
    const bool print_warnings)
{

    namespace interp = ::kynema_sgf::interp;
    const int npts = static_cast<int>(view.pos.size());
    const int npts_vel = static_cast<int>(view.vel_pos.size());
    const std::string who = label.empty() ? "FLLC: " : "FLLC (" + label + "): ";

    if (npts < 2) {
        std::ostringstream msg;
        msg << who << "the correction needs at least 2 actuator points, and "
            << "the line has " << npts;
        amrex::Abort(msg.str());
    }
    if (!(eps_chord > 0.0_rt)) {
        std::ostringstream msg;
        msg << who
            << "epsilon_chord must be larger than 0. Value: " << eps_chord;
        amrex::Abort(msg.str());
    }

    // The velocity points and the force points are two sets of points along
    // the same line, starting at the same place (external turbines: the
    // nodes of the external solver and the actuator points). Each set is
    // located by the distance from its own first point, because the two
    // sets can be positions at different instants of a moving blade.
    data.span_distance_force.resize(npts);
    data.span_distance_vel.resize(npts_vel);
    for (int i = 0; i < npts; ++i) {
        data.span_distance_force[i] = vs::mag(view.pos[i] - view.pos[0]);
    }
    for (int i = 0; i < npts_vel; ++i) {
        data.span_distance_vel[i] = vs::mag(view.vel_pos[i] - view.vel_pos[0]);
    }
    // The same number of points does not make them the same points: compare
    // the distances
    const auto tol = 1.0e-3_rt * data.span_distance_force.back() /
                     static_cast<amrex::Real>(npts);
    data.mismatched_points = npts != npts_vel;
    for (int i = 0; (i < npts) && !data.mismatched_points; ++i) {
        data.mismatched_points =
            std::abs(data.span_distance_vel[i] - data.span_distance_force[i]) >
            tol;
    }
    data.r.resize(npts);
    data.dr.resize(npts);
    data.vel_rel.resize(npts);
    data.optimal_epsilon.resize(npts);
    data.les_velocity.assign(npts, vs::Vector::zero());
    data.optimal_velocity.assign(npts, vs::Vector::zero());
    data.correction_velocity.assign(npts, vs::Vector::zero());
    data.lift.assign(npts, vs::Vector::zero());
    data.grad_lift.assign(npts, vs::Vector::zero());

    for (int i = 0; i < npts; ++i) {
        data.r[i] = vs::mag(view.pos[i] - view.pos[0]);
        data.optimal_epsilon[i] = view.chord[i] * eps_chord;
        // The kernel of the correction divides by the optimal epsilon
        if (!(data.optimal_epsilon[i] > 0.0_rt)) {
            std::ostringstream msg;
            msg << who << "the chord must be larger than 0 at every actuator "
                << "point. Point " << i << " has a chord of " << view.chord[i];
            amrex::Abort(msg.str());
        }
        // The points are located by the distance from the first one
        if ((i > 0) && !(data.r[i] > data.r[i - 1])) {
            std::ostringstream msg;
            msg << who << "the actuator points must move away from the first "
                << "point of the line. Points " << i - 1 << " and " << i
                << " are at distances of " << data.r[i - 1] << " and "
                << data.r[i];
            amrex::Abort(msg.str());
        }
    }

    // The lift is known at the actuator points only, so these points have
    // to resolve the optimal epsilon. Compare it with the distance to the
    // neighbouring points.
    const amrex::Real min_eps_spacing_ratio = 1.0_rt;
    int num_coarse = 0;
    amrex::Real ratio_lowest = std::numeric_limits<amrex::Real>::max();
    for (int i = 0; i < npts; ++i) {
        const auto before = (i > 0) ? data.r[i] - data.r[i - 1] : 0.0_rt;
        const auto after = (i < npts - 1) ? data.r[i + 1] - data.r[i] : 0.0_rt;
        const auto ratio =
            data.optimal_epsilon[i] / amrex::max<amrex::Real>(before, after);
        ratio_lowest = amrex::min<amrex::Real>(ratio_lowest, ratio);
        if (ratio < min_eps_spacing_ratio) {
            ++num_coarse;
        }
    }
    if (print_warnings && (num_coarse > 0)) {
        amrex::OutStream()
            << "\nWARNING: " << who << "the optimal epsilon (chord times "
            << "epsilon_chord) is smaller than the spacing of the actuator "
            << "points at " << num_coarse << " of " << npts
            << " points; the lowest ratio of the two is " << ratio_lowest
            << ". "
            << (data.nonuniform
                    ? "The correction then depends on the number of "
                      "actuator points: use more points or a larger "
                      "epsilon_chord."
                    : "The sum over the actuator points does not resolve "
                      "the correction there: set fllc_nonuniform = true, or "
                      "use more points or a larger epsilon_chord.")
            << std::endl;
    }

    // The correction is handed over at the velocity points. Where these are
    // not the actuator points they have to resolve the optimal epsilon too:
    // the correction can vary over that length, as it does at the tip of a
    // blade, and coarser points show only a few values of it.
    if (data.mismatched_points) {
        int num_coarse_vel = 0;
        amrex::Real ratio_lowest_vel = std::numeric_limits<amrex::Real>::max();
        for (int i = 0; i < npts_vel; ++i) {
            const auto s = data.span_distance_vel[i];
            const auto before =
                (i > 0) ? s - data.span_distance_vel[i - 1] : 0.0_rt;
            const auto after =
                (i < npts_vel - 1) ? data.span_distance_vel[i + 1] - s : 0.0_rt;
            const auto spacing = amrex::max<amrex::Real>(before, after);
            const auto eps_opt =
                interp::linear(data.r, data.optimal_epsilon, s);
            if (eps_opt < min_eps_spacing_ratio * spacing) {
                ++num_coarse_vel;
                ratio_lowest_vel = amrex::min<amrex::Real>(
                    ratio_lowest_vel, eps_opt / spacing);
            }
        }
        if (print_warnings && (num_coarse_vel > 0)) {
            amrex::OutStream()
                << "\nWARNING: " << who << "the optimal epsilon (chord times "
                << "epsilon_chord) is smaller than the spacing of the "
                << "velocity points at " << num_coarse_vel << " of " << npts_vel
                << " points; the lowest ratio of the two is "
                << ratio_lowest_vel << ". The correction can vary over the "
                << "length of the optimal epsilon, as it does at the tip of "
                << "a blade: the velocity handed over at these points then "
                << "depends on where they are. Use more velocity points "
                << "there (the blade nodes of an external solver) or a "
                << "larger epsilon_chord." << std::endl;
        }
    }

    // Central difference
    for (int i = 1; i < npts - 1; ++i) {
        data.dr[i] = vs::mag(view.pos[i + 1] - view.pos[i - 1]) / 2.0_rt;
    }
    data.dr[0] = vs::mag(view.pos[1] - view.pos[0]) / 2.0_rt;
    data.dr[npts - 1] =
        vs::mag(view.pos[npts - 1] - view.pos[npts - 2]) / 2.0_rt;

    //
    //  The following code is used to create a non-uniform distribution of
    //  points. This should follow a spacing given by epsilon/dr and values
    //  above 1 are recommended for converged solutions.
    //
    if (data.nonuniform) {
        // Non-uniform radial distribution variables
        amrex::Real dr;
        amrex::Real dr1;
        amrex::Real eps1;
        amrex::Real dr2;
        amrex::Real eps2;
        // Far more points than any blade needs (a few hundred are typical)
        const amrex::Long max_nonuniform_points = 1000000;

        // Initialize the first location
        data.nonuniform_r.push_back(data.r[0]);
        while (data.nonuniform_r.back() < data.r.back()) {

            amrex::Real r_ =
                data.nonuniform_r.back(); // The latest radial position

            // Interpolate the value of epsilon to the current radial location
            eps1 = interp::linear(data.r, data.optimal_epsilon, r_);
            dr1 = eps1 / data.eps_dr; // spacing needed to match condition

            // Interpolate the value of epsilon to the next radial location
            eps2 = interp::linear(data.r, data.optimal_epsilon, r_ + dr1);
            dr2 = eps2 / data.eps_dr;

            // This will ensure that the spacing will always meet the
            // requirement eps/dr
            dr = amrex::min<amrex::Real>(dr1, dr2);

            // A spacing that does not move the point, or that needs an
            // endless number of points to reach the end of the line
            if (!(r_ + dr > r_) ||
                (data.nonuniform_r.size() >= max_nonuniform_points)) {
                std::ostringstream msg;
                msg << who << "cannot place the non-uniform points: after "
                    << data.nonuniform_r.size() << " points the spacing "
                    << "(chord times epsilon_chord over "
                    << "fllc_epsilon_dr_ratio) is " << dr
                    << " at a distance of " << r_ << " along a line of "
                    << data.r.back() << ". Check the chord, epsilon_chord "
                    << "and fllc_epsilon_dr_ratio.";
                amrex::Abort(msg.str());
            }

            // Append value to the array
            // Ensure that the value is smaller than the tip
            data.nonuniform_r.push_back(
                amrex::min<amrex::Real>(r_ + dr, data.r.back()));
        }

        int npts_n = static_cast<int>(data.nonuniform_r.size());
        data.nonuniform_dr.resize(npts_n);
        data.nonuniform_vel_rel.resize(npts_n);
        data.nonuniform_optimal_epsilon.resize(npts_n);
        data.nonuniform_lift.resize(npts_n);

        // Central difference for non-uniform points
        for (int i = 1; i < npts_n - 1; ++i) {
            data.nonuniform_dr[i] =
                std::abs(data.nonuniform_r[i + 1] - data.nonuniform_r[i - 1]) /
                2.0_rt;
        }
        data.nonuniform_dr[0] =
            std::abs(data.nonuniform_r[1] - data.nonuniform_r[0]) / 2.0_rt;
        data.nonuniform_dr[npts_n - 1] =
            std::abs(
                data.nonuniform_r[npts_n - 1] - data.nonuniform_r[npts_n - 2]) /
            2.0_rt;
    }

    data.initialized = true;
}

void fllc_parse(const utils::ActParser& pp, FLLCData& data)
{
    pp.query("epsilon", data.epsilon);
    pp.query("fllc_relaxation_factor", data.relaxation_factor);
    pp.query("fllc_start_time", data.fllc_start_time);
    std::string typeString = "variable_chord";
    pp.query("fllc_type", typeString);
    const auto match = FLLCTypeMap.count(typeString);
    if (match == 0) {
        amrex::Abort(
            "fllc_parse: Invalid fllc_type provided. Available options are "
            "constant_chord and variable_chord.");
    }
    data.correction_type = FLLCTypeMap.at(typeString);
    pp.query("fllc_nonuniform", data.nonuniform);
    pp.query("fllc_epsilon_dr_ratio", data.eps_dr);

    // With a factor of 0 the correction stays at 0, and above 1 every
    // update overshoots
    if (!(data.relaxation_factor > 0.0_rt) ||
        (data.relaxation_factor > 1.0_rt)) {
        std::ostringstream msg;
        msg << "fllc_parse: fllc_relaxation_factor must be larger than 0 and "
            << "not larger than 1. Value read: " << data.relaxation_factor;
        amrex::Abort(msg.str());
    }
    // The non-uniform points are spaced by the optimal epsilon over this
    // ratio
    if (!(data.eps_dr > 0.0_rt)) {
        std::ostringstream msg;
        msg << "fllc_parse: fllc_epsilon_dr_ratio must be larger than 0. "
            << "Value read: " << data.eps_dr;
        amrex::Abort(msg.str());
    }

    if (!pp.contains("epsilon") || !pp.contains("epsilon_chord")) {
        amrex::Abort(
            "Actuators using the filtered lifting line correction (FLLC) "
            "require specification 'epsilon' and 'epsilon_chord'");
    }
}

} // namespace kynema_sgf::actuator
