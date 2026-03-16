#ifndef RAYLEIGHTAYLORFIELDINIT_H
#define RAYLEIGHTAYLORFIELDINIT_H

#include "AMReX_Array.H"
#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_Geometry.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "AMReX_Gpu.H"

using namespace amrex::literals;

namespace amr_wind {

/** Initialize RayleighTaylor fields
 */
class RayleighTaylorFieldInit
{
    static_assert(
        AMREX_SPACEDIM == 3, "RayleighTaylor requires 3 dimensional mesh");

public:
    RayleighTaylorFieldInit();

    void operator()(
        const amrex::Box& vbx,
        const amrex::Geometry& geom,
        const amrex::Array4<amrex::Real>& density) const;

private:
    //! density value for lower layer
    amrex::Real m_rho_lo{0.5_rt};
    //! density value for upper layer
    amrex::Real m_rho_hi{2.0_rt};
};

} // namespace amr_wind

#endif /* RayleighTaylorFIELDINIT_H */
