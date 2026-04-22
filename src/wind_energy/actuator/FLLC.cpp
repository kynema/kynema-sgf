#include "amr-wind/wind_energy/actuator/FLLC.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator {

void fllc_init(
    FLLCData& data, const ComponentView& view, const amrex::Real eps_chord)
{

    namespace interp = ::amr_wind::interp;
    const int npts = static_cast<int>(view.pos.size());
    data.different_sizes = view.pos.size() != view.vel_pos.size();
    if (data.different_sizes) {
        data.span_distance_force.resize(view.pos.size());
        data.span_distance_vel.resize(view.vel_pos.size());
        for (size_t i = 0; i < view.pos.size(); ++i) {
            data.span_distance_force[i] = vs::mag(view.pos[i]);
        }
        for (size_t i = 0; i < view.vel_pos.size(); ++i) {
            data.span_distance_vel[i] = vs::mag(view.vel_pos[i]);
        }
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
    data.force_point_velocity.assign(npts, vs::Vector::zero());

    for (int i = 0; i < npts; ++i) {
        data.r[i] = vs::mag(view.pos[i] - view.pos[0]);
        data.optimal_epsilon[i] = view.chord[i] * eps_chord;
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
            AMREX_ALWAYS_ASSERT(
                dr > std::numeric_limits<amrex::Real>::epsilon());

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

    if (!pp.contains("epsilon") || !pp.contains("epsilon_chord")) {
        amrex::Abort(
            "Actuators using the filtered lifting line correction (FLLC) "
            "require specification 'epsilon' and 'epsilon_chord'");
    }
}

} // namespace amr_wind::actuator
