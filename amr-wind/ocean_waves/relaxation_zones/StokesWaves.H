#ifndef STOKESWAVES_H
#define STOKESWAVES_H

#include "amr-wind/ocean_waves/relaxation_zones/RelaxationZones.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::ocean_waves {

struct StokesWavesData : public RelaxZonesBaseData
{
    amrex::Real wave_height{0.1_rt};
    amrex::Real wave_length{1.0_rt};
    amrex::Real wave_phase_offset{0.0_rt};
    int order{2};
};

struct StokesWaves : public RelaxZonesType
{
    using InfoType = OceanWavesInfo;
    using MetaType = StokesWavesData;
    using DataType = OceanWavesDataHolder<StokesWaves>;

    static std::string identifier() { return "StokesWaves"; }
};

} // namespace amr_wind::ocean_waves

#endif /* STOKESWAVES_H */
