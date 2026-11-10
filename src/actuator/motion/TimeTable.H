#ifndef ACTUATOR_MOTION_TIME_TABLE_H
#define ACTUATOR_MOTION_TIME_TABLE_H

#include "src/wind_energy/actuator/actuator_types.H"

#include <string>

namespace kynema_sgf::actuator::motion {

/** Piecewise-linear, multi-component time history.
 *
 * Values, derivatives, and integrals use the same interpolation intervals so
 * that motion reconstructed from a velocity history remains consistent.
 */
class TimeTable
{
public:
    void read(
        const std::string& filename,
        int num_values,
        const std::string& extrapolation = "hold");

    [[nodiscard]] bool empty() const { return m_time.empty(); }
    [[nodiscard]] int num_values() const { return m_num_values; }
    [[nodiscard]] const RealList& times() const { return m_time; }
    [[nodiscard]] RealList row(int index) const;

    [[nodiscard]] RealList value(amrex::Real time) const;
    [[nodiscard]] RealList derivative(amrex::Real time) const;
    [[nodiscard]] RealList integral(amrex::Real time) const;

private:
    [[nodiscard]] int interval(amrex::Real time) const;
    void warn_hold(amrex::Real time) const;

    std::string m_filename;
    int m_num_values{0};
    RealList m_time;
    //! Row-major values used directly by interp::linear's ncomp interface
    RealList m_values;
    //! Exact trapezoidal integral of each completed interpolation interval
    amrex::Vector<RealList> m_prefix_integral;
    mutable bool m_warned_below{false};
    mutable bool m_warned_above{false};
    bool m_error_on_extrapolation{false};
};

} // namespace kynema_sgf::actuator::motion

#endif
