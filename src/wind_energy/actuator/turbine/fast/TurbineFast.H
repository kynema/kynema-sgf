#ifndef TURBINEFAST_H
#define TURBINEFAST_H

#include "amr-wind/wind_energy/actuator/turbine/external/TurbineExternal.H"
#include "amr-wind/wind_energy/actuator/turbine/fast/fast_types.H"

namespace amr_wind::actuator {

using TurbineFast =
    TurbineFromExtSolver<::ext_turb::FastTurbine, ::ext_turb::FastSolverData>;

} // namespace amr_wind::actuator

#endif /* TURBINEFAST_H */
