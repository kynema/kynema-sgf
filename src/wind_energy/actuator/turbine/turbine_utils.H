#ifndef TURBINE_UTILS_H
#define TURBINE_UTILS_H

#include "amr-wind/wind_energy/actuator/turbine/turbine_types.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"

namespace amr_wind::actuator::utils {
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

} // namespace amr_wind::actuator::utils

#endif /* TURBINE_UTILS_H */
