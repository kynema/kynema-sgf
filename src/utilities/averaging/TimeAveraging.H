#ifndef TIMEAVERAGING_H
#define TIMEAVERAGING_H

#include "amr-wind/core/Factory.H"
#include "amr-wind/utilities/PostProcessing.H"
#include "amr-wind/utilities/constants.H"

#include "AMReX_Vector.H"

#include <limits>
#include <memory>
#include <string>
#include <map>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

class CFDSim;
class SimTime;
class Field;

namespace averaging {

/** Abstract class for time-averaging of CFD fields.
 *
 *  \ingroup utilities
 */
class FieldTimeAverage
    : public Factory<
          FieldTimeAverage,
          CFDSim&,
          const std::string&,
          const std::string&>
{
public:
    static std::string base_identifier() { return "FieldTimeAverage"; }

    ~FieldTimeAverage() override = default;

    /** Update field averaging at a given timestep
     *
     *  \param time SimTime instance
     *  \param filter_width Time-averaging window specified by user
     *  \param avg_time_interval Time interval at which averaging is performed
     *  \param elapsed_time Time elapsed since averaging was initiated
     */
    virtual void operator()(
        const SimTime& time,
        amrex::Real filter_width,
        amrex::Real avg_time_interval,
        amrex::Real elapsed_time) = 0;

    virtual const std::string& average_field_name() = 0;
};

/** A collection of time-averaged quantities
 */
class TimeAveraging : public PostProcessBase::Register<TimeAveraging>
{
public:
    static std::string identifier() { return "TimeAveraging"; }

    TimeAveraging(CFDSim& /*sim*/, std::string /*label*/);

    ~TimeAveraging() override;

    //! Perform actions before mesh is created
    void pre_init_actions() override;

    void initialize() override;

    void post_advance_work() override;

    void output_actions() override {}

    void post_regrid_actions() override {}

    const std::string& add_averaging(
        const std::string& field_name,
        const std::string& avg_type = "ReAveraging");

private:
    CFDSim& m_sim;

    const std::string m_label;

    //! Vector containing all registered field averagers
    amrex::Vector<std::unique_ptr<FieldTimeAverage>> m_averages;

    //! Fields registered so far to avoid duplication
    std::map<std::string, FieldTimeAverage*> m_registered;

    //! Time to start averaging the fields
    amrex::Real m_start_time{0.0_rt};

    //! Time to stop averaging the fields
    amrex::Real m_stop_time{std::numeric_limits<amrex::Real>::max()};

    //! Time averaging window (in seconds)
    amrex::Real m_filter{amr_wind::constants::LARGE_NUM};

    //! Time averaging time interval (in seconds)
    amrex::Real m_time_interval{-1.};

    //! Accumulated averaging time interval
    amrex::Real m_accumulated_avg_time_interval{0.};
};

} // namespace averaging
} // namespace amr_wind

#endif /* TIMEAVERAGING_H */
