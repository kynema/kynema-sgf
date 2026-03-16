#ifndef RELAX_ZONE_H
#define RELAX_ZONE_H

#include "amr-wind/ocean_waves/OceanWavesTypes.H"
#include "amr-wind/core/vs/vector_space.H"

#include <string>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::ocean_waves {

/** Base data representation of a regular waves
 *
 *  \ingroup ocean_waves
 */
struct RelaxZonesBaseData
{
    amrex::Real zsl{0.0_rt};

    amrex::Real water_depth{0.5_rt};

    // Gravitational constant - wave types that use gravity store value here
    amrex::Real g;

    // Wave generation/absorption parameters
    amrex::Real gen_length{4.0_rt};

    amrex::Real beach_length{8.0_rt};
    amrex::Real beach_length_factor{1.0_rt};

    amrex::Real zone_length_y{0.0_rt};

    amrex::Real current{0.0_rt};

    bool init_wave_field{false};

    bool has_ramp{false};

    bool has_beach{true};

    amrex::Real ramp_period;

    bool regrid_occurred{false};
};

struct RelaxZonesType : public OceanWavesType
{};

} // namespace amr_wind::ocean_waves

#endif // REGULAR_WAVES_H
