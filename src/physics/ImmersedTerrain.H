#ifndef IMMERSEDTERRAIN_H
#define IMMERSEDTERRAIN_H

#include "src/core/Physics.H"
#include "src/core/Field.H"
#include "src/CFDSim.H"

namespace kynema_sgf::immersedterrain {

/** Immersed terrain geometry for the drag-based immersed boundary method
 *  \ingroup physics
 *
 *  Successor to TerrainDrag with a partial (0..1) terrain volume fraction
 *  instead of a binary blanking, so that slopes are not represented as a
 *  staircase. The class only builds geometry fields; the momentum and
 *  temperature forcing consume them.
 *
 *  Fields declared (all with one ghost cell):
 *
 *  - ``terrain_fraction`` (Real, 1 component): fraction of the cell volume
 *    occupied by terrain. 0 in fluid, 1 inside terrain, in between for cells
 *    cut by the terrain surface: the fraction of the cell column below the
 *    terrain height at the cell center,
 *    \f[ \beta = \min\left(1, \max\left(0,
 *         \frac{h - z_{bottom}}{\Delta z}\right)\right) \f]
 *
 *  - ``terrain_mask`` (Int, 1 component): cell classification
 *      - 0 ``mask_fluid``   : fluid cell away from the terrain
 *      - 1 ``mask_solid``   : cell fully inside the terrain
 *      - 2 ``mask_surface`` : a partially filled cell, or a fluid cell sharing
 *        a face with a wall cell on any of its six sides. The six-way search
 *        catches the vertical faces of steep terrain and buildings that a
 *        purely vertical search misses.
 *    Use with ``FieldRefinement``: ``field_error = 1.5`` tags only the
 *    surface band, ``field_error = 0.5`` tags surface and solid cells.
 *
 *  - ``terrain_surface`` (Real, 1 component): terrain height h(x,y) at the
 *    cell center, used for the wall distance of the no-slip flux.
 *
 *  - ``terrain_drag_rate`` (Real, 1 component), only with
 *    ``implicit_projection = true``: the immersed drag rate of each cell,
 *    \f$ C = w C_d / \Delta z \f$ with the drag weight \f$ w \f$ of
 *    immersed_wall::solid_weight. Its presence makes the nodal and MAC
 *    projections apply the drag together with the pressure,
 *    \f[ \mathbf{u}^{n+1} = \frac{\mathbf{u}^{**} - \Delta t \nabla p / \rho}
 *         {1 + C \Delta t} \quad \Rightarrow \quad
 *         \nabla \cdot \left( \frac{\Delta t}{\rho (1 + C \Delta t)}
 *         \nabla p \right) = \nabla \cdot \frac{\mathbf{u}^{**}}
 *         {1 + C \Delta t} \f]
 *    so the terrain acts as a fluid of density \f$ \rho (1 + C \Delta t) \f$
 *    and the pressure gradient produces no velocity inside it. The same
 *    effective density is the time coefficient of the implicit diffusion
 *    solve, which keeps momentum from diffusing a depth
 *    \f$ \sqrt{\nu \Delta t} \f$ into the body every step.
 *    ImmersedDragForcing then skips its explicit drag. Without this option
 *    the velocity inside the body is first order in dt.
 *
 *  - ``terrain_diffusion_xf/yf/zf`` (Real, face centered): factor applied to
 *    the face diffusion coefficients of every equation solved with the shared
 *    diffusion operator. Wall faces are scaled by dx_f / d1 so that the
 *    discrete flux mu (u - 0)/dx_f equals the flux to a no-slip wall at the
 *    true wall position,
 *    \f[ q = \mu \frac{u}{d_1}, \qquad d_1 = \left\{ \begin{array}{ll}
 *         z_c - h & \mbox{center} \\
 *         (1 - \beta) \Delta z / 2 & \mbox{fraction, partial cell}
 *         \end{array} \right. \f]
 *    With ``fraction`` weighting the partial cell stands for the centroid of
 *    its fluid part, so the face above it spans
 *    \f$ (2 - \beta) \Delta z / 2 \f$ and is scaled by
 *    \f$ 2 / (2 - \beta) \f$. Side walls use d1 = dx_f / 2.
 */
class ImmersedTerrain : public Physics::Register<ImmersedTerrain>
{
public:
    static std::string identifier() { return "ImmersedTerrain"; }

    //! Values stored in the ``terrain_mask`` int field
    static constexpr int mask_fluid = 0;
    static constexpr int mask_solid = 1;
    static constexpr int mask_surface = 2;

    //! Component indices of the ``terrain_surface`` field
    static constexpr int surf_height = 0;

    explicit ImmersedTerrain(CFDSim& sim);

    ~ImmersedTerrain() override = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void pre_init_actions() override {}

    void post_init_actions() override {}

    void post_regrid_actions() override;

    void pre_advance_work() override {}

    void post_advance_work() override {}

private:
    CFDSim& m_sim;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;

    //! Terrain volume fraction per cell (0..1)
    Field& m_terrain_fraction;
    //! Cell classification: fluid / solid / surface
    IntField& m_terrain_mask;
    //! Terrain height at the cell center
    Field& m_terrain_surface;

    //! Terrain file (x, y, h on a flat grid)
    std::string m_terrain_file{"terrain.amrwind"};

    //! With center weighting, a cell is solid at or above this fraction
    amrex::Real m_solid_threshold{0.5_rt};
    //! ``fraction`` (drag weight beta) or ``center`` (1 when the cell center
    //! is inside the terrain, else 0); shared with the source terms
    std::string m_drag_weight{"fraction"};

    //! Apply the immersed drag implicitly through the projections
    bool m_implicit_projection{false};

    //! turbulence.model = Laminar; with fraction weighting the partial cell
    //! is then a fluid cell and only solid cells take the drag
    bool m_laminar{false};
    //! Drag coefficient C_d shared with ImmersedDragForcing
    amrex::Real m_drag_coefficient{10.0_rt};
    //! Drag rate w C_d / dz, declared only when m_implicit_projection is true
    Field* m_terrain_drag_rate{nullptr};

    //! Face factors for the diffusion coefficients (no-slip wall flux)
    amrex::Array<Field*, AMREX_SPACEDIM> m_diffusion_factor{
        {nullptr, nullptr, nullptr}};
};

} // namespace kynema_sgf::immersedterrain

#endif
