#include "src/wind_energy/actuator/wing/FlatPlate.H"
#include "src/wind_energy/actuator/wing/flat_plate_ops.H"
#include "src/wind_energy/actuator/ActuatorModel.H"

namespace kynema_sgf::actuator {

template class ActModel<FlatPlate, ActSrcLine>;

} // namespace kynema_sgf::actuator
