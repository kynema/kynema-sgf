#include "src/actuator/disk/UniformCt.H"
#include "src/actuator/disk/Joukowsky.H"
#include "src/actuator/disk/uniform_ct_ops.H"
#include "src/actuator/disk/Joukowsky_ops.H"
#include "src/actuator/ActuatorModel.H"
#include "src/actuator/disk/disk_spreading.H"

namespace kynema_sgf::actuator {
template class ActModel<UniformCt, ActSrcDisk>;
template class ActModel<Joukowsky, ActSrcDisk>;

} // namespace kynema_sgf::actuator
