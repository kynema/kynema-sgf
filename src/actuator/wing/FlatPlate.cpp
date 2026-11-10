#include "src/actuator/wing/FlatPlate.H"
#include "src/actuator/wing/flat_plate_ops.H"
#include "src/actuator/ActuatorModel.H"

namespace kynema_sgf::actuator {

template class ActModel<FlatPlate, ActSrcLine>;

} // namespace kynema_sgf::actuator
