#ifndef LINEARWAVES_H
#define LINEARWAVES_H

#include "src/ocean_waves/relaxation_zones/RelaxationZones.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::ocean_waves {

struct LinearWavesData : public RelaxZonesBaseData
{
    amrex::Real wave_height{0.1_rt};
    amrex::Real wave_length{1.0_rt};
    amrex::Real wave_phase_offset{0.0_rt};
    // Remove the mass transport of the waves with a uniform return current
    // (Fenton's u2-bar = 0, default); mean_current is the resulting Eulerian
    // mean current (Fenton's u1-bar) applied to the wave kinematics. When
    // false, the kinematics are those for zero Eulerian mean current.
    bool zero_mass_transport{true};
    amrex::Real mean_current{0.0_rt};
};

struct LinearWaves : public RelaxZonesType
{
    using InfoType = OceanWavesInfo;
    using MetaType = LinearWavesData;
    using DataType = OceanWavesDataHolder<LinearWaves>;

    static std::string identifier() { return "LinearWaves"; }
};

} // namespace kynema_sgf::ocean_waves

#endif /* LINEARWAVES_H */
