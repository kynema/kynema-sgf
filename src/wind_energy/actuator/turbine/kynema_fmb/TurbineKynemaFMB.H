#ifndef TURBINEKYNEMAFMB_H
#define TURBINEKYNEMAFMB_H

#include "src/wind_energy/actuator/turbine/external/TurbineExternal.H"
#include "src/wind_energy/actuator/turbine/kynema_fmb/kynema_fmb_types.H"

namespace kynema_sgf::actuator {

using TurbineKynemaFMB = TurbineFromExtSolver<
    ::ext_turb::KynemaFMBTurbine,
    ::ext_turb::KynemaFMBSolverData>;

} // namespace kynema_sgf::actuator

#endif /* TURBINEKYNEMAFMB_H */
