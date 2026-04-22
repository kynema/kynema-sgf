#ifndef FIELDDESCTYPES_H
#define FIELDDESCTYPES_H

#include <cstdint>

namespace amr_wind {

/** Possible time states for a given Field instance
 *  \ingroup fields
 *
 *  All defined fields have at least one state, `New`, as their default state.
 *  For state variables, there are currently two states `New` and `Old`. `Old` .
 *  indicates the latest solution and `New` indicates the solution at the next
 *  timestep the solution will advance to.
 */
enum class FieldState : std::uint8_t {
    NP1 = 0,   ///< Latest state, also aliased to `New`
    N = 1,     ///< Previous state, also aliased to `Old`
    NM1 = 2,   ///< State at `(n - 1)` th timestep
    NPH = 3,   ///< State at `(n + 1/2)` (intermediate) timestep
    NMH = 4,   ///< State at `(n - 1/2)` timestep
    New = NP1, ///< Same as FieldState::NP1
    Old = N,   ///< Same as FieldState::N
};

/** Possible locations for a field in a structured mesh
 *  \ingroup fields
 */
enum class FieldLoc : std::uint8_t {
    CELL,  ///< Cell-centered (default)
    NODE,  ///< Node-centered (e.g., for pressure)
    XFACE, ///< Face-centered in x-direction (e.g., face normal velocity)
    YFACE, ///< Face-centered in y-direction
    ZFACE  ///< Face-centered in z-direction
};

//! Coarse-to-fine field interpolation options
//! \ingroup fields
enum class FieldInterpolator : std::uint8_t {
    PiecewiseConstant, ///< Constant across cell
    CellConsLinear,    ///< Linear interpolation
    NodeBilinear,      ///< Bilinear nodal interpolation
    FaceDivFree,       ///< Divergence free face interpolation
    FaceLinear         ///< Linear face interpolation
};

} // namespace amr_wind

#endif /* FIELDDESCTYPES_H */
