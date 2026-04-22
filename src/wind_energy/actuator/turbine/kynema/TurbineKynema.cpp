#include "amr-wind/wind_energy/actuator/turbine/kynema/TurbineKynema.H"
#include "amr-wind/wind_energy/actuator/turbine/kynema/turbine_kynema_ops.H"
#include "amr-wind/wind_energy/actuator/ActuatorModel.H"

namespace amr_wind::actuator {

template class ActModel<TurbineKynema, ActSrcLine>;
template class ActModel<TurbineKynema, ActSrcDisk>;

} // namespace amr_wind::actuator

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
