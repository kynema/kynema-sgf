#include "amr-wind/wind_energy/actuator/disk/UniformCt.H"
#include "amr-wind/wind_energy/actuator/disk/Joukowsky.H"
#include "amr-wind/wind_energy/actuator/disk/uniform_ct_ops.H"
#include "amr-wind/wind_energy/actuator/disk/Joukowsky_ops.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"
#include "amr-wind/wind_energy/actuator/disk/disk_spreading.H"

namespace amr_wind::actuator {
template class ActModel<UniformCt, ActSrcDisk>;
template class ActModel<Joukowsky, ActSrcDisk>;

} // namespace amr_wind::actuator
