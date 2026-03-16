#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "amr-wind/immersed_boundary/IBTypes.H"
#include "amr-wind/core/vs/vector_space.H"

#include <string>
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::ib {

/** Base data representation of an IB geometry.
 *
 *  \ingroup immersed_boundary
 */
struct BluffBodyBaseData
{
    bool is_moving{false};

    bool is_mms{false};

    bool has_wall_model{false};

    amrex::Vector<amrex::Real> vel_bc{0.0_rt, 0.0_rt, 0.0_rt};

    //! Total integrated forces on the immersed body
    amrex::Vector<amrex::Real> frc{0.0_rt, 0.0_rt, 0.0_rt};
};

struct BluffBodyType : public IBType
{};

} // namespace amr_wind::ib

#endif /* GEOMETRY_H */
