#ifndef PARTIALTERRAINDRAG_H
#define PARTIALTERRAINDRAG_H

#include "src/core/Physics.H"
#include "src/core/Field.H"
#include "src/CFDSim.H"

namespace kynema_sgf::partialterraindrag {

/** Partial Terrain Drag physics with volume fraction blanking
 *  \ingroup physics
 *
 *  This module implements terrain blanking using volume fractions (0.0 to 1.0)
 *  instead of binary blanking (0 or 1). This provides smoother transitions
 *  at terrain boundaries and more accurate representation of terrain geometry.
 */

class PartialTerrainDrag : public Physics::Register<PartialTerrainDrag>
{
public:
    static std::string identifier() { return "PartialTerrainDrag"; }

    explicit PartialTerrainDrag(CFDSim& sim);

    ~PartialTerrainDrag() override = default;

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

    //! Partial blanking field for terrain (volume fraction: 0.0 to 1.0)
    Field& m_terrain_blank;
    //! Terrain drag force term (volume fraction: 0.0 to 1.0)
    Field& m_terrain_drag;

    //! Terrain file
    std::string m_terrain_file{"terrain.amrwind"};

    //! Roughness file
    std::string m_roughness_file{"terrain.roughness"};

    //! Roughness fields
    Field& m_terrainz0;
    Field& m_terrain_height;
    Field& m_terrain_damping;

    //! Lateral Damping - X (East)
    amrex::Real m_damp_east_slope{0.0_rt};
    amrex::Real m_damp_east_full{0.0_rt};
    amrex::Real m_horizontal_tau{20.0_rt};
    amrex::Real m_horizontal_abl_height{1e15_rt};
    amrex::Real m_horizontal_slope_end{1e30_rt};
    amrex::Real m_vertical_slope{0.0_rt};
    amrex::Real m_vertical_full{1e30_rt};

    //! Lateral Damping - West (-X)
    amrex::Real m_damp_west_slope{0.0_rt};
    amrex::Real m_damp_west_full{0.0_rt};

    //! Lateral Damping - South (-Y)
    amrex::Real m_damp_south_slope{0.0_rt};
    amrex::Real m_damp_south_full{0.0_rt};

    //! Lateral Damping - North (Y)
    amrex::Real m_damp_north_slope{0.0_rt};
    amrex::Real m_damp_north_full{0.0_rt};

    //! Blanking method selection
    //! Options: "volume_fraction" (default) or "distance_function"
    std::string m_blanking_method{"volume_fraction"};

    //! Smoothing length for distance function approach (in grid cells)
    //! Controls the transition width at terrain boundary
    amrex::Real m_smoothing_length{1.0_rt};
};

} // namespace kynema_sgf::partialterraindrag

#endif
