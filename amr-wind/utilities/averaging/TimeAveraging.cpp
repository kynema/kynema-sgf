#include <utility>

#include "amr-wind/utilities/averaging/TimeAveraging.H"
#include "amr-wind/utilities/averaging/ReAveraging.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/CFDSim.H"

#include "AMReX_ParmParse.H"

namespace amr_wind::averaging {

TimeAveraging::TimeAveraging(CFDSim& sim, std::string label)
    : m_sim(sim), m_label(std::move(label))
{}

TimeAveraging::~TimeAveraging() = default;

void TimeAveraging::pre_init_actions()
{
    //! Different averaging types
    amrex::Vector<std::string> labels;
    {
        amrex::ParmParse pp(m_label);
        pp.getarr("labels", labels);
        ioutils::assert_with_message(
            ioutils::all_distinct(labels),
            "Duplicates in " + m_label + ".labels");
        pp.query("averaging_start_time", m_start_time);
        pp.query("averaging_stop_time", m_stop_time);
        pp.query("averaging_time_interval", m_time_interval);
        pp.query("averaging_window", m_filter);
    }

    amrex::Print() << "TimeAveraging: Initializing " << m_label << '\n';

    for (const auto& lbl : labels) {
        //! Fields to be averaged
        amrex::Vector<std::string> fnames;
        std::string avg_type;
        const std::string pp_key = m_label + "." + lbl;
        amrex::ParmParse pp1(pp_key);
        pp1.getarr("fields", fnames);
        ioutils::assert_with_message(
            ioutils::all_distinct(fnames),
            "Duplicates in " + pp_key + ".fields");
        pp1.get("averaging_type", avg_type);

        amrex::Print() << "    - initializing average labeled " << lbl
                       << ", type " << avg_type << '\n';

        for (const auto& fname : fnames) {
            const std::string key = fname + "_" + avg_type;

            // Create the averaging entity
            m_averages.emplace_back(
                FieldTimeAverage::create(avg_type, m_sim, m_label, fname));

            // Track fields that have an average
            m_registered.emplace(key, m_averages.back().get());
        }
    }
}

void TimeAveraging::initialize() {}

const std::string& TimeAveraging::add_averaging(
    const std::string& field_name, const std::string& avg_type)
{
    const std::string key = field_name + "_" + avg_type;
    const auto found = m_registered.find(key);

    // If the field was already registered then just return the field
    if (found != m_registered.end()) {
        return m_registered[key]->average_field_name();
    }

    // Create and register new average
    m_averages.emplace_back(
        FieldTimeAverage::create(avg_type, m_sim, m_label, field_name));
    return m_averages.back()->average_field_name();
}

void TimeAveraging::post_advance_work()
{
    const auto& time = m_sim.time();
    const auto cur_time = time.new_time();
    const auto cur_dt = time.delta_t();

    m_accumulated_avg_time_interval += cur_dt;

    // Check the following:
    //   1. if we are phase averaging (i.e., only accumulating the average at a
    //   certain frequency), check to see if we are on a time step in which
    //   averaging should be performed.  This is useful for simulations with
    //   periodic behavior, such as regular waves or actuator lines rotating at
    //   fixed rotor speed.
    //   2. if we are within the averaging time period requested by the user
    const auto t_tol = amr_wind::constants::LOOSE_TOL * cur_dt;
    const bool do_phase_avg =
        (m_time_interval < 0. ||
         ((cur_time - m_start_time + t_tol) / m_time_interval) -
                 std::floor(
                     (cur_time - m_start_time + t_tol) / m_time_interval) <
             cur_dt / m_time_interval);
    const bool do_avg =
        (((cur_time >= m_start_time) && (cur_time < m_stop_time)) &&
         do_phase_avg);
    if (!do_avg) {
        return;
    }

    const amrex::Real elapsed_time = (cur_time - m_start_time);
    for (const auto& avg : m_averages) {
        (*avg)(time, m_filter, m_accumulated_avg_time_interval, elapsed_time);
    }
    m_accumulated_avg_time_interval = 0.;
}

} // namespace amr_wind::averaging
