#ifndef TURBINEFAST_H
#define TURBINEFAST_H

#include "src/wind_energy/actuator/turbine/external/TurbineExternal.H"
#include "src/wind_energy/actuator/turbine/fast/fast_types.H"

namespace kynema_sgf::actuator {

using TurbineFast =
    TurbineFromExtSolver<::ext_turb::FastTurbine, ::ext_turb::FastSolverData>;

} // namespace kynema_sgf::actuator

#endif /* TURBINEFAST_H */
