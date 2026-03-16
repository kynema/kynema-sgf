#include "amr-wind/physics/TerrainDrag.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "amr-wind/utilities/trig_ops.H"
#include "amr-wind/utilities/IOManager.H"
#include "amr-wind/utilities/io_utils.H"
#include "amr-wind/utilities/linear_interpolation.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::terraindrag {

namespace {} // namespace

TerrainDrag::TerrainDrag(CFDSim& sim)
    : m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_terrain_blank(sim.repo().declare_int_field("terrain_blank", 1, 1, 1))
    , m_terrain_drag(sim.repo().declare_int_field("terrain_drag", 1, 1, 1))
    , m_terrainz0(sim.repo().declare_field("terrainz0", 1, 1, 1))
    , m_terrain_height(sim.repo().declare_field("terrain_height", 1, 1, 1))
{

    m_terrain_is_waves = sim.physics_manager().contains("OceanWaves") &&
                         !sim.repo().field_exists("vof");

    if (!m_terrain_is_waves) {
        amrex::ParmParse pp(identifier());
        pp.query("terrain_file", m_terrain_file);
        pp.query("roughness_file", m_roughness_file);
    } else {
        m_wave_volume_fraction = &m_repo.get_field(m_wave_volume_fraction_name);
        m_wave_negative_elevation =
            &m_repo.get_field(m_wave_negative_elevation_name);
    }

    m_sim.io_manager().register_output_int_var("terrain_drag");
    m_sim.io_manager().register_output_int_var("terrain_blank");
    m_sim.io_manager().register_io_var("terrainz0");
    m_sim.io_manager().register_io_var("terrain_height");

    m_terrain_blank.setVal(0);
    m_terrain_drag.setVal(0);
    m_terrainz0.set_default_fillpatch_bc(m_sim.time());
    m_terrain_height.set_default_fillpatch_bc(m_sim.time());
}

void TerrainDrag::initialize_fields(int level, const amrex::Geometry& geom)
{
    if (m_terrain_is_waves) {
        return;
    }

    BL_PROFILE("amr-wind::" + this->identifier() + "::initialize_fields");

    //! Reading the Terrain Coordinates from  file
    amrex::Vector<amrex::Real> xterrain;
    amrex::Vector<amrex::Real> yterrain;
    amrex::Vector<amrex::Real> zterrain;
    ioutils::read_flat_grid_file(m_terrain_file, xterrain, yterrain, zterrain);

    // No checks for the file as it is optional currently
    amrex::Vector<amrex::Real> xrough;
    amrex::Vector<amrex::Real> yrough;
    amrex::Vector<amrex::Real> z0rough;
    std::ifstream file(m_roughness_file, std::ios::in);
    if (file.good()) {
        ioutils::read_flat_grid_file(m_roughness_file, xrough, yrough, z0rough);
    }
    file.close();

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    auto& blanking = m_terrain_blank(level);
    auto& terrainz0 = m_terrainz0(level);
    auto& terrain_height = m_terrain_height(level);
    auto& drag = m_terrain_drag(level);
    const auto xterrain_size = xterrain.size();
    const auto yterrain_size = yterrain.size();
    const auto zterrain_size = zterrain.size();
    amrex::Gpu::DeviceVector<amrex::Real> d_xterrain(xterrain_size);
    amrex::Gpu::DeviceVector<amrex::Real> d_yterrain(yterrain_size);
    amrex::Gpu::DeviceVector<amrex::Real> d_zterrain(zterrain_size);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, xterrain.begin(), xterrain.end(),
        d_xterrain.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, yterrain.begin(), yterrain.end(),
        d_yterrain.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, zterrain.begin(), zterrain.end(),
        d_zterrain.begin());
    const auto* xterrain_ptr = d_xterrain.data();
    const auto* yterrain_ptr = d_yterrain.data();
    const auto* zterrain_ptr = d_zterrain.data();
    // Copy Roughness to gpu
    const auto xrough_size = xrough.size();
    const auto yrough_size = yrough.size();
    const auto z0rough_size = z0rough.size();
    amrex::Gpu::DeviceVector<amrex::Real> d_xrough(xrough_size);
    amrex::Gpu::DeviceVector<amrex::Real> d_yrough(yrough_size);
    amrex::Gpu::DeviceVector<amrex::Real> d_z0rough(z0rough_size);
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, xrough.begin(), xrough.end(),
        d_xrough.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, yrough.begin(), yrough.end(),
        d_yrough.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, z0rough.begin(), z0rough.end(),
        d_z0rough.begin());
    const auto* xrough_ptr = d_xrough.data();
    const auto* yrough_ptr = d_yrough.data();
    const auto* z0rough_ptr = d_z0rough.data();
    auto levelBlanking = blanking.arrays();
    auto levelDrag = drag.arrays();
    auto levelz0 = terrainz0.arrays();
    auto levelheight = terrain_height.arrays();
    amrex::ParallelFor(
        blanking, m_terrain_blank.num_grow(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real x = prob_lo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = prob_lo[1] + ((j + 0.5_rt) * dx[1]);
            const amrex::Real z = prob_lo[2] + ((k + 0.5_rt) * dx[2]);
            const amrex::Real terrainHt = interp::bilinear(
                xterrain_ptr, xterrain_ptr + xterrain_size, yterrain_ptr,
                yterrain_ptr + yterrain_size, zterrain_ptr, x, y);
            levelBlanking[nbx](i, j, k, 0) =
                static_cast<int>((z <= terrainHt) && (z > prob_lo[2]));
            levelheight[nbx](i, j, k, 0) = terrainHt;

            amrex::Real roughz0 = 0.1_rt;
            if (xrough_size > 0) {
                roughz0 = interp::bilinear(
                    xrough_ptr, xrough_ptr + xrough_size, yrough_ptr,
                    yrough_ptr + yrough_size, z0rough_ptr, x, y);
            }
            levelz0[nbx](i, j, k, 0) = roughz0;
        });
    amrex::Gpu::streamSynchronize();
    amrex::ParallelFor(
        blanking, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            if ((levelBlanking[nbx](i, j, k, 0) == 0) && (k > 0) &&
                (levelBlanking[nbx](i, j, k - 1, 0) == 1)) {
                levelDrag[nbx](i, j, k, 0) = 1;
            } else {
                levelDrag[nbx](i, j, k, 0) = 0;
            }
        });
    amrex::Gpu::streamSynchronize();
}

void TerrainDrag::post_init_actions()
{
    if (!m_terrain_is_waves) {
        return;
    }
    BL_PROFILE("amr-wind::" + this->identifier() + "::post_init_actions");
    convert_waves_to_terrain_fields();
}

void TerrainDrag::pre_advance_work()
{
    if (!m_terrain_is_waves) {
        return;
    }
    BL_PROFILE("amr-wind::" + this->identifier() + "::pre_advance_work");
    convert_waves_to_terrain_fields();
}

void TerrainDrag::post_regrid_actions()
{
    if (m_terrain_is_waves) {
        convert_waves_to_terrain_fields();
    } else {
        const int nlevels = m_sim.repo().num_active_levels();
        for (int lev = 0; lev < nlevels; ++lev) {
            initialize_fields(lev, m_sim.repo().mesh().Geom(lev));
        }
    }
}

void TerrainDrag::convert_waves_to_terrain_fields()
{
    const int nlevels = m_sim.repo().num_active_levels();
    // Uniform, low roughness for waves
    m_terrainz0.setVal(1.0e-4_rt);
    for (int level = 0; level < nlevels; ++level) {
        const auto geom = m_sim.repo().mesh().Geom(level);
        const auto& dx = geom.CellSizeArray();
        const auto& prob_lo = geom.ProbLoArray();
        auto& blanking = m_terrain_blank(level);
        auto& terrain_height = m_terrain_height(level);
        auto& drag = m_terrain_drag(level);

        auto levelBlanking = blanking.arrays();
        auto levelDrag = drag.arrays();
        auto levelHeight = terrain_height.arrays();

        const auto negative_wave_elevation =
            (*m_wave_negative_elevation)(level).const_arrays();
        const auto wave_vol_frac =
            (*m_wave_volume_fraction)(level).const_arrays();

        // Get terrain blanking from ocean waves fields
        amrex::ParallelFor(
            blanking, m_terrain_blank.num_grow(),
            [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                const amrex::Real z = prob_lo[2] + ((k + 0.5_rt) * dx[2]);
                levelBlanking[nbx](i, j, k, 0) = static_cast<int>(
                    (wave_vol_frac[nbx](i, j, k) >= 0.5_rt) &&
                    (z > prob_lo[2]));
                levelHeight[nbx](i, j, k, 0) =
                    -negative_wave_elevation[nbx](i, j, k);
            });
        amrex::ParallelFor(
            blanking, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
                if ((levelBlanking[nbx](i, j, k, 0) == 0) && (k > 0) &&
                    (levelBlanking[nbx](i, j, k - 1, 0) == 1)) {
                    levelDrag[nbx](i, j, k, 0) = 1;
                } else {
                    levelDrag[nbx](i, j, k, 0) = 0;
                }
            });
    }
}

} // namespace amr_wind::terraindrag
