#include "amr-wind/wind_energy/actuator/wing/FixedWing.H"
#include "amr-wind/wind_energy/actuator/wing/fixed_wing_ops.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"

namespace amr_wind::actuator {

template class ActModel<FixedWing, ActSrcLine>;

} // namespace amr_wind::actuator
