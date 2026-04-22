#ifndef TURBINEKYNEMA_H
#define TURBINEKYNEMA_H

#include "amr-wind/wind_energy/actuator/turbine/external/TurbineExternal.H"
#include "amr-wind/wind_energy/actuator/turbine/kynema/kynema_types.H"

namespace amr_wind::actuator {

using TurbineKynema = TurbineFromExtSolver<
    ::ext_turb::KynemaTurbine,
    ::ext_turb::KynemaSolverData>;

} // namespace amr_wind::actuator

#endif /* TURBINEKYNEMA_H */
