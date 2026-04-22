#ifndef INCFLO_ENUMS_H
#define INCFLO_ENUMS_H

#include <cstdint>

enum struct DiffusionType : std::uint8_t {
    Invalid,
    Explicit,
    Crank_Nicolson,
    Implicit
};

enum struct BC : std::uint8_t {
    pressure_outflow,
    mass_inflow,
    mass_inflow_outflow,
    zero_gradient,
    no_slip_wall,
    slip_wall,
    wall_model,
    fixed_gradient,
    periodic,
    wave_generation,
    symmetric_wall,
    undefined
};

#endif /* INCFLO_ENUMS_H */
