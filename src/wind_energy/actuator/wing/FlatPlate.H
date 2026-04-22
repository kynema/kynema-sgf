#ifndef FLATPLATE_H
#define FLATPLATE_H

#include "amr-wind/wind_energy/actuator/wing/ActuatorWing.H"
#include "amr-wind/wind_energy/actuator/aero/AirfoilTable.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::actuator {

struct FlatPlateData : public WingBaseData
{
    amrex::Real chord_inp{1.0_rt};
    amrex::Real cd_factor{0.0_rt};

    ThinAirfoil aflookup;
};

struct FlatPlate : public WingType
{
    using InfoType = ActInfo;
    using GridType = ActGrid;
    using MetaType = FlatPlateData;
    using DataType = ActDataHolder<FlatPlate>;

    static std::string identifier() { return "FlatPlate"; }
};

} // namespace amr_wind::actuator

#endif /* FLATPLATE_H */
