#include "src/wind_energy/actuator/sector/ActuatorSector.H"
#include "src/wind_energy/actuator/sector/actuator_sector_ops.H"
#include "src/wind_energy/actuator/ActuatorModel.H"

namespace kynema_sgf::actuator {

template class ActModel<ActuatorSector, ActSrcSector>;

} // namespace kynema_sgf::actuator
