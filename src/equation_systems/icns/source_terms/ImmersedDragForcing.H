#ifndef IMMERSEDDRAGFORCING_H
#define IMMERSEDDRAGFORCING_H

#include "src/equation_systems/icns/MomentumSource.H"
#include "src/core/SimTime.H"
#include "src/CFDSim.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::pde::icns {

/** Momentum forcing for terrain represented by the ImmersedTerrain physics
 *
 *  \ingroup icns_src
 *
 *  Consumes ``terrain_fraction`` declared by ImmersedTerrain. Immersed drag
 *  in every cell with a non-zero drag weight w, relaxing the velocity toward
 *  zero,
 *  \f[ \frac{\partial \mathbf{u}}{\partial t} = \dots - C \mathbf{u},
 *      \qquad C = \frac{w C_d}{\Delta z} \f]
 *  The relaxation is integrated exactly over the time step: the source is
 *  \f$ -C_{eff} \mathbf{u}^n \f$ with
 *  \f[ C_{eff} = \frac{1 - e^{-C \Delta t}}{\Delta t} \quad \Rightarrow
 *      \quad \mathbf{u}^{n+1} = \mathbf{u}^n e^{-C \Delta t} \f]
 *  so it is stable and monotone for any C and dt with no ad hoc limiters
 *  (\f$ C_{eff} \to C \f$ for small \f$ C \Delta t \f$ and
 *  \f$ C_{eff} \to 1 / \Delta t \f$ for large).
 *
 *  The drag weight follows ``ImmersedTerrain.drag_weight``: ``fraction`` uses
 *  w = beta, ``center`` uses w = 1 where the cell center is inside the
 *  terrain. With laminar flow and ``fraction`` weighting the partial cell is
 *  a fluid cell whose wall is carried by the no-slip wall flux, so only cells
 *  entirely inside the terrain take the drag.
 *
 *  With ``ImmersedTerrain.implicit_projection`` the drag is applied inside the
 *  projections and the diffusion solve instead and this source adds nothing.
 *  As an explicit source the drag is not seen by the implicit diffusion solve,
 *  which leaves a splitting error of depth \f$ \sqrt{\nu \Delta t} \f$ at the
 *  wall: first order, and converging only when dt is refined with the mesh.
 */
class ImmersedDragForcing : public MomentumSource::Register<ImmersedDragForcing>
{
public:
    static std::string identifier() { return "ImmersedDragForcing"; }

    explicit ImmersedDragForcing(const CFDSim& sim);

    ~ImmersedDragForcing() override;

    void operator()(
        int lev, FieldState fstate, amrex::MultiFab& src_term) const override;

private:
    const SimTime& m_time;
    const CFDSim& m_sim;
    const amrex::AmrCore& m_mesh;
    const Field& m_velocity;

    //! Immersed drag relaxation rate is drag_coefficient / dz
    amrex::Real m_drag_coefficient{10.0_rt};
    //! ``fraction`` or ``center`` drag weight, read from ImmersedTerrain
    std::string m_drag_weight{"fraction"};
    //! Terrain fraction at or above which ``center`` treats a cell as solid
    amrex::Real m_solid_threshold{0.5_rt};
    //! turbulence.model = Laminar
    bool m_is_laminar{false};
};

} // namespace kynema_sgf::pde::icns

#endif
