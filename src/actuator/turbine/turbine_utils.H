#ifndef TURBINE_UTILS_H
#define TURBINE_UTILS_H

#include <string>
#include "src/actuator/turbine/turbine_types.H"
#include "src/actuator/ActParser.H"

namespace kynema_sgf::actuator::utils {

/* Warning for FLLC inputs that differ by direction, empty when none do
 *
 *  The FLLC takes the first value of epsilon and epsilon_min and the second
 *  of epsilon_chord (inputs in the order chord, thickness, span).
 */
std::string fllc_epsilon_warning(
    const std::string& label,
    const vs::Vector& eps_inp,
    const vs::Vector& eps_min,
    const vs::Vector& eps_chord);

void read_inputs(
    TurbineBaseData& tdata, TurbineInfo& info, const utils::ActParser& pp);

void prepare_netcdf_file(
    const std::string& /*ncfile*/,
    const TurbineBaseData& /*meta*/,
    const TurbineInfo& /*info*/,
    const ActGrid& /*grid*/);

void write_netcdf(
    const std::string& /*ncfile*/,
    const TurbineBaseData& /*meta*/,
    const TurbineInfo& /*info*/,
    const ActGrid& /*grid*/,
    amrex::Real /*time*/);

} // namespace kynema_sgf::actuator::utils

#endif /* TURBINE_UTILS_H */
