#ifndef TERRAINDRAG_H
#define TERRAINDRAG_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"

namespace amr_wind::terraindrag {

/** Terraindrag Flow physics
 *  \ingroup physics
 */

class TerrainDrag : public Physics::Register<TerrainDrag>
{
public:
    static std::string identifier() { return "TerrainDrag"; }

    explicit TerrainDrag(CFDSim& sim);

    ~TerrainDrag() override = default;

    void initialize_fields(int level, const amrex::Geometry& geom) override;

    void pre_init_actions() override {}

    void post_init_actions() override;

    void post_regrid_actions() override;

    void pre_advance_work() override;

    void post_advance_work() override {}

    //! Terrain Drag for waves

    void convert_waves_to_terrain_fields();

    [[nodiscard]] std::string wave_velocity_field_name() const
    {
        return m_wave_velocity_name;
    }

    [[nodiscard]] std::string wave_negative_elevation_name() const
    {
        return m_wave_negative_elevation_name;
    }

private:
    CFDSim& m_sim;
    const FieldRepo& m_repo;
    const amrex::AmrCore& m_mesh;

    //! Blanking field for terrain or buildings
    IntField& m_terrain_blank;
    //! Terrain drag force term
    IntField& m_terrain_drag;

    //! Terrain file
    std::string m_terrain_file{"terrain.amrwind"};

    //! Roughness file
    std::string m_roughness_file{"terrain.roughness"};

    //! Roughness fields
    Field& m_terrainz0;
    Field& m_terrain_height;
    Field& m_terrain_damping;

    //! Lateral Damping - X
    amrex::Real m_damp_east_slope{0};
    amrex::Real m_damp_east_full{0};
    amrex::Real m_horizontal_tau{20.0};
    amrex::Real m_horizontal_abl_height{1e15};
    amrex::Real m_horizontal_slope_end{1e30};
    amrex::Real m_horizontal_free_atmosphere{1e30};
    amrex::Real m_vertical_slope{0};
    amrex::Real m_vertical_full{1e30};
    //! Lateral Damping - (-X)
    amrex::Real m_damp_west_slope{0};
    amrex::Real m_damp_west_full{0};
    //! Lateral Damping - (-Y)
    amrex::Real m_damp_south_slope{0};
    amrex::Real m_damp_south_full{0};
    //! Lateral Damping - Y
    amrex::Real m_damp_north_slope{0};
    amrex::Real m_damp_north_full{0};

    //! Terrain Drag for waves
    bool m_terrain_is_waves{false};
    Field* m_wave_volume_fraction{nullptr};
    Field* m_wave_negative_elevation{nullptr};
    const std::string m_wave_volume_fraction_name{"ow_vof"};
    const std::string m_wave_negative_elevation_name{"ow_levelset"};
    const std::string m_wave_velocity_name{"ow_velocity"};
};
} // namespace amr_wind::terraindrag

#endif
