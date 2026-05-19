#include "src/wind_energy/actuator/turbine/kynema_fmb/TurbineKynemaFMB.H"
#include "src/wind_energy/actuator/turbine/kynema_fmb/turbine_kynema_fmb_ops.H"
#include "src/wind_energy/actuator/ActuatorModel.H"

namespace kynema_sgf::actuator {

template class ActModel<TurbineKynemaFMB, ActSrcLine>;
template class ActModel<TurbineKynemaFMB, ActSrcDisk>;

} // namespace kynema_sgf::actuator

namespace ext_turb {
template <>
std::string ext_id<KynemaFMBTurbine>()
{
    return "TurbineKynemaFMB";
}
template <>
std::string ext_id<KynemaFMBSolverData>()
{
    return "KynemaFMB";
}
} // namespace ext_turb
