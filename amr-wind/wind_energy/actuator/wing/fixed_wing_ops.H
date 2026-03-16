#ifndef FIXED_WING_OPS_H
#define FIXED_WING_OPS_H

#include <algorithm>

#include "amr-wind/wind_energy/actuator/wing/FixedWing.H"
#include "amr-wind/wind_energy/actuator/actuator_ops.H"
#include "amr-wind/wind_energy/actuator/wing/wing_ops.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator::ops {

template <>
struct ReadInputsOp<FixedWing, ActSrcLine>
{
    void operator()(FixedWing::DataType& data, const utils::ActParser& pp)
    {
        auto& wdata = data.meta();
        auto& info = data.info();
        pp.get("num_points", wdata.num_pts);
        pp.get("start", wdata.start);
        pp.get("end", wdata.end);
        pp.query_either("epsilon", wdata.eps_inp);
        pp.query_either("epsilon_chord", wdata.epsilon_chord);
        if (!pp.contains("pitch_timetable")) {
            // Constant pitch is required without timetable
            pp.get("pitch", wdata.pitch);
        }
        pp.get("airfoil_table", wdata.airfoil_file);
        pp.query("airfoil_type", wdata.airfoil_type);
        pp.queryarr("span_locs", wdata.span_locs);
        pp.queryarr("chord", wdata.chord_inp);
        bool use_fllc = false;
        pp.query("fllc", use_fllc);

        // If spanwise components of Gaussian should be ignored
        pp.query("disable_spanwise_gaussian", wdata.gauss_2D);
        // If spanwise components of Gaussian and forces should be normalized
        // (on by default, not used unless gauss_2D = true)
        pp.query("normalize_spanwise", wdata.normalize_2D_spanwise);
        // If certain components of force should be ignored
        // (check that correct number of components are specified)
        amrex::Vector<int> force_coord_flags_query_inp;
        pp.queryarr("active_force_dirs", force_coord_flags_query_inp);
        if (!force_coord_flags_query_inp.empty()) {
            AMREX_ALWAYS_ASSERT(
                force_coord_flags_query_inp.size() == AMREX_SPACEDIM);
            wdata.force_coord_flags = force_coord_flags_query_inp;
        }

        // If velocity magnitude for force calculation should be prescribed by
        // user, not measured from the flow
        pp.query("prescribed_uinf", wdata.prescribed_uinf);

        // Initialize tables for pitch actuation (in degrees)
        pp.query("pitch_timetable", wdata.pitch_timetable_file);
        if (!wdata.pitch_timetable_file.empty()) {
            std::ifstream ifh(wdata.pitch_timetable_file, std::ios::in);
            if (!ifh.good()) {
                amrex::Abort(
                    "Cannot find input file: " + wdata.pitch_timetable_file);
            }
            amrex::Real data_time;
            amrex::Real data_pitch_deg;
            ifh.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            while (ifh >> data_time) {
                ifh >> data_pitch_deg;
                wdata.time_table.push_back(data_time);
                wdata.pitch_table.push_back(data_pitch_deg);
            }
        }

        pp.query("motion_type", wdata.motion_type);
        {
            // do nothing for default, "none"
            if (amrex::toLower(wdata.motion_type) == "linear") {
                // get velocity
                pp.get("velocity", wdata.vel_tr);
            } else if (amrex::toLower(wdata.motion_type) == "sine") {
                // get parameters
                pp.get("sine_period", wdata.s_period);
                pp.get("sine_vector", wdata.s_vector);
            } else if (amrex::toLower(wdata.motion_type) != "none") {
                amrex::Abort(
                    "Invalid FixedWingLine motion type. Valid options are "
                    "none, linear, and sine.");
            }
        }
        if (use_fllc) {
            wdata.fllc = std::make_unique<FLLCData>();
            fllc_parse(pp, *(wdata.fllc));
        }

        if (!pp.contains("epsilon") && !pp.contains("epsilon_chord")) {
            amrex::Abort(
                "Actuator fixed wing requires specification of one or both "
                "of 'epsilon' or 'epsilon_chord'");
        }

        AMREX_ALWAYS_ASSERT(wdata.span_locs.size() == wdata.chord_inp.size());

        amrex::Real max_chord = *std::ranges::max_element(wdata.chord_inp);
        amrex::Real max_eps = *std::ranges::max_element(wdata.eps_inp);
        amrex::Real max_epsc = *std::ranges::max_element(wdata.epsilon_chord);
        amrex::Real search_radius =
            amrex::max(max_eps, max_epsc) * max_chord * 3.0_rt;
        const auto& p1 = wdata.start;
        const auto& p2 = wdata.end;
        // Bounding box quantities used for all motion types
        amrex::Real minpx = amrex::min(p1.x(), p2.x());
        amrex::Real maxpx = amrex::max(p1.x(), p2.x());
        amrex::Real minpy = amrex::min(p1.y(), p2.y());
        amrex::Real maxpy = amrex::max(p1.y(), p2.y());
        amrex::Real minpz = amrex::min(p1.z(), p2.z());
        amrex::Real maxpz = amrex::max(p1.z(), p2.z());
        // Bounding box limits for 2D actuator
        if (wdata.gauss_2D) {
            // Assume largest component of span needs entire domain, but the
            // others do not need the bounding box to grow
            const auto wspan = wdata.end - wdata.start;
            // Flag for checking periodicity in spanwise direction
            bool warn_per{false};
            std::string sdir;
            auto per = data.sim().mesh().Geom()[0].periodicity();
            if (std::abs(wspan.x()) >=
                amrex::max<amrex::Real>(
                    std::abs(wspan.y()), std::abs(wspan.z()))) {
                minpx = data.sim().mesh().Geom(0).ProbLoArray()[0];
                maxpx = data.sim().mesh().Geom(0).ProbHiArray()[0];
                warn_per = !per.isPeriodic(0);
                sdir = "x";
            } else if (std::abs(wspan.y()) >= std::abs(wspan.z())) {
                minpy = data.sim().mesh().Geom(0).ProbLoArray()[1];
                maxpy = data.sim().mesh().Geom(0).ProbHiArray()[1];
                warn_per = !per.isPeriodic(1);
                sdir = "y";
            } else {
                minpz = data.sim().mesh().Geom(0).ProbLoArray()[2];
                maxpz = data.sim().mesh().Geom(0).ProbHiArray()[2];
                warn_per = !per.isPeriodic(2);
                sdir = "z";
            }
            if (warn_per) {
                amrex::Print()
                    << "\nWARNING: fixed_wing_ops: Detected spanwise direction "
                    << sdir
                    << " is not periodic, though 2D Gaussian is being "
                       "used.\n\n";
            }
        }
        // Set up bounding box depending on motion type
        if (amrex::toLower(wdata.motion_type) == "none") {
            info.bound_box = amrex::RealBox(
                minpx - search_radius, minpy - search_radius,
                minpz - search_radius, maxpx + search_radius,
                maxpy + search_radius, maxpz + search_radius);
        } else if (amrex::toLower(wdata.motion_type) == "linear") {
            // Extend bounding box in case of velocity
            constexpr amrex::Real tiny = std::numeric_limits<float>::epsilon();
            // If velocity is present, assume wing can travel entire domain
            const bool up_ = wdata.vel_tr.x() > tiny;
            const amrex::Real up_ext =
                data.sim().mesh().Geom(0).ProbHiArray()[0];
            const bool um_ = wdata.vel_tr.x() < -tiny;
            const amrex::Real um_ext =
                data.sim().mesh().Geom(0).ProbLoArray()[0];
            const bool vp_ = wdata.vel_tr.y() > tiny;
            const amrex::Real vp_ext =
                data.sim().mesh().Geom(0).ProbHiArray()[1];
            const bool vm_ = wdata.vel_tr.y() < -tiny;
            const amrex::Real vm_ext =
                data.sim().mesh().Geom(0).ProbLoArray()[1];
            const bool wp_ = wdata.vel_tr.z() > tiny;
            const amrex::Real wp_ext =
                data.sim().mesh().Geom(0).ProbHiArray()[2];
            const bool wm_ = wdata.vel_tr.z() < -tiny;
            const amrex::Real wm_ext =
                data.sim().mesh().Geom(0).ProbLoArray()[2];
            info.bound_box = amrex::RealBox(
                um_ ? um_ext : minpx - search_radius,
                vm_ ? vm_ext : minpy - search_radius,
                wm_ ? wm_ext : minpz - search_radius,
                up_ ? up_ext : maxpx + search_radius,
                vp_ ? vp_ext : maxpy + search_radius,
                wp_ ? wp_ext : maxpz + search_radius);
        } else if (amrex::toLower(wdata.motion_type) == "sine") {
            info.bound_box = amrex::RealBox(
                minpx - wdata.s_vector.x() - search_radius,
                minpy - wdata.s_vector.y() - search_radius,
                minpz - wdata.s_vector.z() - search_radius,
                maxpx + wdata.s_vector.x() + search_radius,
                maxpy + wdata.s_vector.y() + search_radius,
                maxpz + wdata.s_vector.z() + search_radius);
        }
    }
};

template <>
struct InitDataOp<FixedWing, ActSrcLine>
{
    void operator()(FixedWing::DataType& data)
    {
        wing::init_data_structures(data.meta(), data.grid());

        auto& meta = data.meta();
        {
            const int npts = data.meta().num_pts;
            auto& grid = data.grid();
            // Process chord information
            const auto wlen = vs::mag(grid.pos.back() - grid.pos.front());
            RealList wx(npts);
            for (int i = 0; i < npts; ++i) {
                wx[i] = vs::mag(grid.pos[i] - grid.pos[0]) / wlen;
            }
            meta.chord.resize(npts);
            ::amr_wind::interp::linear_monotonic(
                meta.span_locs, meta.chord_inp, wx, meta.chord);
            meta.epsilon_chord = {
                meta.epsilon_chord.x(), meta.epsilon_chord.z(),
                meta.epsilon_chord.y()};
            meta.eps_inp = {
                meta.eps_inp.x(), meta.eps_inp.z(), meta.eps_inp.y()};
            for (int i = 0; i < npts; ++i) {
                for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                    const auto eps = meta.epsilon_chord[n] * meta.chord[i];
                    grid.epsilon[i][n] =
                        amrex::max<amrex::Real>(meta.eps_inp[n], eps);
                }
            }
            // Copy Gaussian flags to grid struct if 2D
            if (meta.gauss_2D) {
                // Local coords are chord, span, thickness internally
                grid.dcoord_flags = vs::Vector(1.0_rt, 0.0_rt, 1.0_rt);
            }
        }

        meta.aflookup =
            AirfoilLoader::load_airfoil(meta.airfoil_file, meta.airfoil_type);
    }
};

template <>
AMREX_INLINE const AirfoilTable&
airfoil_lookup<FixedWing>(FixedWing::DataType& data)
{
    return *data.meta().aflookup;
}

} // namespace amr_wind::actuator::ops

#endif /* FIXED_WING_OPS_H */
