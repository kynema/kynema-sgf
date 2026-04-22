#ifndef TURBINEKYNEMA_H
#define TURBINEKYNEMA_H

#include "src/wind_energy/actuator/turbine/external/TurbineExternal.H"
#include "src/wind_energy/actuator/turbine/kynema/kynema_types.H"

namespace kynema_sgf::actuator {

using TurbineKynema = TurbineFromExtSolver<
    ::ext_turb::KynemaTurbine,
    ::ext_turb::KynemaSolverData>;

} // namespace kynema_sgf::actuator

#endif /* TURBINEKYNEMA_H */
