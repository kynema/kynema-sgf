#include "src/actuator/turbine/fast/TurbineFast.H"
#include "src/actuator/turbine/fast/turbine_fast_ops.H"
#include "src/actuator/ActuatorModel.H"

namespace kynema_sgf::actuator {

template class ActModel<TurbineFast, ActSrcLine>;
template class ActModel<TurbineFast, ActSrcDisk>;

} // namespace kynema_sgf::actuator

namespace ext_turb {
template <>
std::string ext_id<FastTurbine>()
{
    return "TurbineFast";
}
template <>
std::string ext_id<FastSolverData>()
{
    return "OpenFAST";
}
} // namespace ext_turb
