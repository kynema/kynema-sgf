#ifndef LINEARWAVES_H
#define LINEARWAVES_H

#include "amr-wind/ocean_waves/relaxation_zones/RelaxationZones.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::ocean_waves {

struct LinearWavesData : public RelaxZonesBaseData
{
    amrex::Real wave_height{0.1_rt};
    amrex::Real wave_length{1.0_rt};
    amrex::Real wave_phase_offset{0.0_rt};
};

struct LinearWaves : public RelaxZonesType
{
    using InfoType = OceanWavesInfo;
    using MetaType = LinearWavesData;
    using DataType = OceanWavesDataHolder<LinearWaves>;

    static std::string identifier() { return "LinearWaves"; }
};

} // namespace amr_wind::ocean_waves

#endif /* LINEARWAVES_H */
