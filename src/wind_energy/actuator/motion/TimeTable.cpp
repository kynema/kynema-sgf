#include "src/wind_energy/actuator/motion/TimeTable.H"

#include "src/utilities/linear_interpolation.H"

#include <fstream>
#include <limits>

#include "AMReX.H"
#include "AMReX_Print.H"

using namespace amrex::literals;

namespace kynema_sgf::actuator::motion {

void TimeTable::read(
    const std::string& filename,
    const int num_values,
    const std::string& extrapolation)
{
    m_filename = filename;
    m_num_values = num_values;
    m_time.clear();
    m_values.clear();
    m_prefix_integral.clear();
    m_warned_below = false;
    m_warned_above = false;
    const auto extrapolation_lower = amrex::toLower(extrapolation);
    if ((extrapolation_lower != "hold") && (extrapolation_lower != "error")) {
        amrex::Abort("Actuator timetable extrapolation must be hold or error");
    }
    m_error_on_extrapolation = (extrapolation_lower == "error");
    std::ifstream stream(filename);
    if (!stream.good()) {
        amrex::Abort("Cannot open actuator motion timetable: " + filename);
    }

    stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    amrex::Real time;
    while (stream >> time) {
        RealList row(num_values);
        for (auto& item : row) {
            if (!(stream >> item)) {
                amrex::Abort(
                    "Invalid actuator motion timetable row in: " + filename);
            }
        }
        if (!m_time.empty() && time <= m_time.back()) {
            amrex::Abort(
                "Actuator motion timetable times must be strictly "
                "increasing: " +
                filename);
        }
        m_time.push_back(time);
        m_values.insert(m_values.end(), row.begin(), row.end());
    }
    if (m_time.empty()) {
        amrex::Abort("Actuator motion timetable contains no data: " + filename);
    }

    m_prefix_integral.assign(m_time.size(), RealList(num_values, 0.0_rt));
    for (int i = 1; i < static_cast<int>(m_time.size()); ++i) {
        const amrex::Real dt = m_time[i] - m_time[i - 1];
        for (int n = 0; n < num_values; ++n) {
            m_prefix_integral[i][n] =
                m_prefix_integral[i - 1][n] +
                0.5_rt * dt *
                    (m_values[(i - 1) * m_num_values + n] +
                     m_values[i * m_num_values + n]);
        }
    }
}

RealList TimeTable::row(const int index) const
{
    AMREX_ALWAYS_ASSERT(index >= 0 && index < static_cast<int>(m_time.size()));
    const auto begin = m_values.begin() + index * m_num_values;
    return RealList(begin, begin + m_num_values);
}

void TimeTable::warn_hold(const amrex::Real time) const
{
    if (m_error_on_extrapolation &&
        ((time < m_time.front()) || (time > m_time.back()))) {
        amrex::Abort(
            "Actuator motion timetable evaluated outside its time range: " +
            m_filename);
    }
    if (time < m_time.front() && !m_warned_below) {
        amrex::Print() << "WARNING: Holding first value of actuator motion "
                          "timetable '"
                       << m_filename << "' before time " << m_time.front()
                       << ".\n";
        m_warned_below = true;
    }
    if (time > m_time.back() && !m_warned_above) {
        amrex::Print() << "WARNING: Holding last value of actuator motion "
                          "timetable '"
                       << m_filename << "' after time " << m_time.back()
                       << ".\n";
        m_warned_above = true;
    }
}

int TimeTable::interval(const amrex::Real time) const
{
    if (m_time.size() == 1 || time <= m_time.front()) {
        return 0;
    }
    if (time >= m_time.back()) {
        return static_cast<int>(m_time.size()) - 2;
    }
    return ::kynema_sgf::interp::bisection_search(
               m_time.begin(), m_time.end(), time)
        .idx;
}

RealList TimeTable::value(const amrex::Real time) const
{
    warn_hold(time);
    RealList result(m_num_values);
    for (int n = 0; n < m_num_values; ++n) {
        result[n] = ::kynema_sgf::interp::linear(
            m_time, m_values, time, m_num_values, n);
    }
    return result;
}

RealList TimeTable::derivative(const amrex::Real time) const
{
    warn_hold(time);
    RealList result(m_num_values, 0.0_rt);
    if (m_time.size() == 1 || time < m_time.front() || time > m_time.back()) {
        return result;
    }
    const int i = interval(time);
    const amrex::Real dt = m_time[i + 1] - m_time[i];
    for (int n = 0; n < m_num_values; ++n) {
        result[n] = (m_values[(i + 1) * m_num_values + n] -
                     m_values[i * m_num_values + n]) /
                    dt;
    }
    return result;
}

RealList TimeTable::integral(const amrex::Real time) const
{
    warn_hold(time);
    RealList result(m_num_values, 0.0_rt);
    if (time <= m_time.front()) {
        const amrex::Real dt = time - m_time.front();
        for (int n = 0; n < m_num_values; ++n) {
            result[n] = dt * m_values[n];
        }
        return result;
    }
    if (time >= m_time.back()) {
        result = m_prefix_integral.back();
        const amrex::Real dt = time - m_time.back();
        for (int n = 0; n < m_num_values; ++n) {
            result[n] += dt * m_values[(m_time.size() - 1) * m_num_values + n];
        }
        return result;
    }
    const int i = interval(time);
    result = m_prefix_integral[i];
    const amrex::Real dt = time - m_time[i];
    const amrex::Real interval_dt = m_time[i + 1] - m_time[i];
    for (int n = 0; n < m_num_values; ++n) {
        const amrex::Real slope = (m_values[(i + 1) * m_num_values + n] -
                                   m_values[i * m_num_values + n]) /
                                  interval_dt;
        result[n] +=
            m_values[i * m_num_values + n] * dt + 0.5_rt * slope * dt * dt;
    }
    return result;
}

} // namespace kynema_sgf::actuator::motion
