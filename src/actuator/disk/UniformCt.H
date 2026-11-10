#ifndef UNIFORMCT_H
#define UNIFORMCT_H

#include "src/wind_energy/actuator/disk/ActuatorDisk.H"

namespace kynema_sgf::actuator {

struct UniformCtData : public DiskBaseData
{};

struct UniformCt : public DiskType
{
    using InfoType = ActInfo;
    using GridType = ActGrid;
    using MetaType = UniformCtData;
    using DataType = ActDataHolder<UniformCt>;

    static std::string identifier() { return "UniformCt"; }
};

} // namespace kynema_sgf::actuator

#endif /* UNIFORMCT_H */
