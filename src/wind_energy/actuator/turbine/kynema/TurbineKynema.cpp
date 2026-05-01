#include "src/wind_energy/actuator/turbine/kynema/TurbineKynema.H"
#include "src/wind_energy/actuator/turbine/kynema/turbine_kynema_ops.H"
#include "src/wind_energy/actuator/ActuatorModel.H"

namespace kynema_sgf::actuator {

template class ActModel<TurbineKynema, ActSrcLine>;
template class ActModel<TurbineKynema, ActSrcDisk>;

} // namespace kynema_sgf::actuator

namespace ext_turb {
template <>
std::string ext_id<KynemaTurbine>()
{
    return "TurbineKynema";
}
template <>
std::string ext_id<KynemaSolverData>()
{
    return "Kynema";
}
} // namespace ext_turb
