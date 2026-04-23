#include "src/wind_energy/actuator/wing/FixedWing.H"
#include "src/wind_energy/actuator/wing/fixed_wing_ops.H"
#include "src/wind_energy/actuator/ActuatorModel.H"

namespace kynema_sgf::actuator {

template class ActModel<FixedWing, ActSrcLine>;

} // namespace kynema_sgf::actuator
