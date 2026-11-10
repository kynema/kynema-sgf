#include "src/actuator/wing/FixedWing.H"
#include "src/actuator/wing/fixed_wing_ops.H"
#include "src/actuator/ActuatorModel.H"

namespace kynema_sgf::actuator {

template class ActModel<FixedWing, ActSrcLine>;

} // namespace kynema_sgf::actuator
