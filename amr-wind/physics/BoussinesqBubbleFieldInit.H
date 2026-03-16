#ifndef BOUSSINESQBUBBLEFIELDINIT_H
#define BOUSSINESQBUBBLEFIELDINIT_H

#include "AMReX_Array.H"
#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_Geometry.H"
#include "AMReX_REAL.H"
#include "AMReX_Vector.H"
#include "AMReX_Gpu.H"

using namespace amrex::literals;

namespace amr_wind {

/** Initialize BoussinesqBubble fields
 */
class BoussinesqBubbleFieldInit
{
    static_assert(
        AMREX_SPACEDIM == 3, "BoussinesqBubble requires 3 dimensional mesh");

public:
    BoussinesqBubbleFieldInit();

    void operator()(
        const amrex::Box& vbx,
        const amrex::Geometry& geom,
        const amrex::Array4<amrex::Real>& velocity,
        const amrex::Array4<amrex::Real>& density,
        const amrex::Array4<amrex::Real>& tracer) const;

private:
    //! Initial bubble location
    amrex::Vector<amrex::Real> m_loc{0.25_rt, 0.25_rt, 0.5_rt};

    //! Initial density field
    amrex::Real m_rho{1.0_rt};

    //! tracer inner value
    amrex::Real m_tracer_inner{0.0_rt};

    //! tracer outer value
    amrex::Real m_tracer_outer{0.01_rt};

    //! tracer bubble radius value
    amrex::Real m_tracer_radius{0.1_rt};
};

} // namespace amr_wind

#endif /* BoussinesqBubbleFIELDINIT_H */
