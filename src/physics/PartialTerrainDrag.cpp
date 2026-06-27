#include "src/physics/PartialTerrainDrag.H"
#include "src/CFDSim.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "src/utilities/trig_ops.H"
#include "src/utilities/IOManager.H"
#include "src/utilities/io_utils.H"
#include "src/utilities/linear_interpolation.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::partialterraindrag {

namespace {} // namespace

PartialTerrainDrag::PartialTerrainDrag(CFDSim& sim)
    : m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_terrain_blank(
          sim.repo().declare_field("partial_terrain_blank", 1, 1, 1))
    , m_terrain_drag(sim.repo().declare_field("partial_terrain_drag", 1, 1, 1))
    , m_terrainz0(sim.repo().declare_field("partial_terrainz0", 1, 1, 1))
    , m_terrain_height(
          sim.repo().declare_field("partial_terrain_height", 1, 1, 1))
    , m_terrain_damping(
          sim.repo().declare_field("partial_terrain_damping", 1, 1, 1))
{
    amrex::ParmParse pp(identifier());
    pp.query("terrain_file", m_terrain_file);
    pp.query("roughness_file", m_roughness_file);

    m_sim.io_manager().register_output_var("partial_terrain_drag");
    m_sim.io_manager().register_output_var("partial_terrain_blank");
    m_sim.io_manager().register_io_var("partial_terrainz0");
    m_sim.io_manager().register_io_var("partial_terrain_height");
    m_sim.io_manager().register_io_var("partial_terrain_damping");

    m_terrain_blank.setVal(0.0);
    m_terrain_drag.setVal(0.0);
    m_terrain_damping.setVal(0.0);
    m_terrainz0.set_default_fillpatch_bc(m_sim.time());
    m_terrain_height.set_default_fillpatch_bc(m_sim.time());
    m_terrain_damping.set_default_fillpatch_bc(m_sim.time());

    pp.query("damp_east_slope", m_damp_east_slope);
    pp.query("damp_east_full", m_damp_east_full);
    pp.query("damp_west_slope", m_damp_west_slope);
    pp.query("damp_west_full", m_damp_west_full);
    pp.query("damp_north_slope", m_damp_north_slope);
    pp.query("damp_north_full", m_damp_north_full);
    pp.query("damp_south_slope", m_damp_south_slope);
    pp.query("damp_south_full", m_damp_south_full);
    pp.query("horizontal_time_scale", m_horizontal_tau);
    pp.query("horizontal_abl_height", m_horizontal_abl_height);
    pp.query("horizontal_slope_end", m_horizontal_slope_end);
    pp.query("vertical_slope", m_vertical_slope);
    pp.query("vertical_full", m_vertical_full);
    pp.query("blanking_method", m_blanking_method);
    pp.query("smoothing_length", m_smoothing_length);
}

void PartialTerrainDrag::initialize_fields(
    int level, const amrex::Geometry& geom)
{
    BL_PROFILE("kynema-sgf::" + this->identifier() + "::initialize_fields");

    //! Reading the Terrain Coordinates from file
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
    const auto& prob_hi = geom.ProbHiArray();
    auto& blanking = m_terrain_blank(level);
    auto& terrainz0 = m_terrainz0(level);
    auto& terrain_height = m_terrain_height(level);
    auto& drag = m_terrain_drag(level);
    auto& damping = m_terrain_damping(level);

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
    auto levelDamping = damping.arrays();

    const auto dlo = geom.Domain().smallEnd();
    const auto dhi = geom.Domain().bigEnd();

    // Determine blanking method
    const bool use_distance_function =
        (m_blanking_method == "distance_function");
    const amrex::Real smooth_len =
        m_smoothing_length * dx[2]; // Convert to physical length

    // Calculate partial blanking using selected method
    amrex::ParallelFor(
        blanking, m_terrain_blank.num_grow(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real x = prob_lo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = prob_lo[1] + ((j + 0.5_rt) * dx[1]);
            const amrex::Real z = prob_lo[2] + ((k + 0.5_rt) * dx[2]);

            const amrex::Real terrainHt = interp::bilinear(
                xterrain_ptr, xterrain_ptr + xterrain_size, yterrain_ptr,
                yterrain_ptr + yterrain_size, zterrain_ptr, x, y);

            // Store terrain height for all cells
            levelheight[nbx](i, j, k, 0) = terrainHt;

            amrex::Real volume_fraction = 0.0_rt;

            if (use_distance_function) {
                // Distance function approach: smooth transition using tanh
                // Distance from cell center to terrain surface
                const amrex::Real dist = z - terrainHt;

                // Smooth blanking function: 0.5 * (1 - tanh(dist / smooth_len))
                // Values: 1.0 below terrain, 0.5 at terrain, 0.0 above terrain
                volume_fraction =
                    0.5_rt * (1.0_rt - std::tanh(dist / smooth_len));

                // Clamp to [0, 1]
                volume_fraction =
                    amrex::max(0.0_rt, amrex::min(1.0_rt, volume_fraction));
            } else {
                // Volume fraction approach: linear interpolation within cell
                // Cell boundaries
                const amrex::Real z_bottom = prob_lo[2] + (k * dx[2]);
                const amrex::Real z_top = prob_lo[2] + ((k + 1) * dx[2]);

                // Fully below terrain
                if (z_top <= terrainHt) {
                    volume_fraction = 1.0_rt;
                }
                // Partially intersecting terrain
                else if (z_bottom < terrainHt && z_top > terrainHt) {
                    // Linear approximation: fraction of cell below terrain
                    volume_fraction = (terrainHt - z_bottom) / dx[2];
                    // Clamp to [0, 1]
                    volume_fraction =
                        amrex::max(0.0_rt, amrex::min(1.0_rt, volume_fraction));
                }
                // Fully above terrain
                else {
                    volume_fraction = 0.0_rt;
                }
            }

            // Only blank cells above ground level
            if (z > prob_lo[2]) {
                levelBlanking[nbx](i, j, k, 0) = volume_fraction;
            } else {
                levelBlanking[nbx](i, j, k, 0) = 0.0_rt;
            }

            // Set roughness
            amrex::Real roughz0 = 0.1_rt;
            if (xrough_size > 0) {
                roughz0 = interp::bilinear(
                    xrough_ptr, xrough_ptr + xrough_size, yrough_ptr,
                    yrough_ptr + yrough_size, z0rough_ptr, x, y);
            }
            levelz0[nbx](i, j, k, 0) = roughz0;
        });
    amrex::Gpu::streamSynchronize();

    // Calculate drag field: cells just above or adjacent to blanked cells get
    // drag force For partial blanking, we use a smooth transition
    amrex::ParallelFor(
        blanking, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real current_blank = levelBlanking[nbx](i, j, k, 0);
            amrex::Real drag_fraction = 0.0_rt;

            amrex::Real max_neighbor_blank = 0.0_rt;
            if (k > dlo[2]) {
                max_neighbor_blank = amrex::max(
                    max_neighbor_blank, levelBlanking[nbx](i, j, k - 1, 0));
            }
            if (k < dhi[2]) {
                max_neighbor_blank = amrex::max(
                    max_neighbor_blank, levelBlanking[nbx](i, j, k + 1, 0));
            }
            if (i > dlo[0]) {
                max_neighbor_blank = amrex::max(
                    max_neighbor_blank, levelBlanking[nbx](i - 1, j, k, 0));
            }
            if (i < dhi[0]) {
                max_neighbor_blank = amrex::max(
                    max_neighbor_blank, levelBlanking[nbx](i + 1, j, k, 0));
            }
            if (j > dlo[1]) {
                max_neighbor_blank = amrex::max(
                    max_neighbor_blank, levelBlanking[nbx](i, j - 1, k, 0));
            }
            if (j < dhi[1]) {
                max_neighbor_blank = amrex::max(
                    max_neighbor_blank, levelBlanking[nbx](i, j + 1, k, 0));
            }

            // If current cell has low blanking and any neighbor has high
            // blanking
            if (current_blank < 0.5_rt && max_neighbor_blank > 0.5_rt) {
                // Full drag at interface
                drag_fraction = 1.0_rt;
            } else if (current_blank > 0.0_rt && current_blank < 1.0_rt) {
                // Partial drag for partially blanked cells
                drag_fraction = 1.0_rt - current_blank;
            }

            levelDrag[nbx](i, j, k, 0) = drag_fraction;
        });
    amrex::Gpu::streamSynchronize();

    // Lateral damping setup (same as original implementation)
    const amrex::Real horizontal_tau = m_horizontal_tau;
    const amrex::Real horizontal_abl_height = m_horizontal_abl_height;
    const amrex::Real z_sloped = m_horizontal_slope_end;
    const amrex::Real vertical_slope = m_vertical_slope;
    const amrex::Real vertical_full = m_vertical_full;
    const amrex::Real damping_east_start =
        prob_hi[0] - (m_damp_east_full + m_damp_east_slope);
    const amrex::Real damping_east_end = prob_hi[0] - m_damp_east_full;
    // West
    const amrex::Real damping_west_start =
        prob_lo[0] + (m_damp_west_full + m_damp_west_slope);
    const amrex::Real damping_west_end = prob_lo[0] + m_damp_west_full;
    // North
    const amrex::Real damping_north_start =
        prob_hi[1] - (m_damp_north_full + m_damp_north_slope);
    const amrex::Real damping_north_end = prob_hi[1] - m_damp_north_full;
    // South
    const amrex::Real damping_south_start =
        prob_lo[1] + (m_damp_south_full + m_damp_south_slope);
    const amrex::Real damping_south_end = prob_lo[1] + m_damp_south_full;

    amrex::ParallelFor(
        damping, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
            amrex::Real horizontal_coeff_east = 0.0_rt;
            amrex::Real horizontal_coeff_north = 0.0_rt;
            amrex::Real horizontal_coeff_west = 0.0_rt;
            amrex::Real horizontal_coeff_south = 0.0_rt;
            amrex::Real vertical_coeff = 0.0_rt;
            const amrex::Real x = prob_lo[0] + (i + 0.5_rt) * dx[0];
            const amrex::Real y = prob_lo[1] + (j + 0.5_rt) * dx[1];
            const amrex::Real z = prob_lo[2] + (k + 0.5_rt) * dx[2];

            // East damping
            if (x < damping_east_start) {
                horizontal_coeff_east = 0.0_rt;
            } else if (x >= damping_east_end) {
                horizontal_coeff_east = 1.0_rt;
            } else {
                const amrex::Real term = std::sin(
                    std::numbers::pi_v<amrex::Real> * 0.5_rt *
                    (x - damping_east_start) /
                    (damping_east_end - damping_east_start));
                horizontal_coeff_east = term * term;
            }

            // West damping
            if (x > damping_west_start) {
                horizontal_coeff_west = 0.0_rt;
            } else if (x <= damping_west_end) {
                horizontal_coeff_west = 1.0_rt;
            } else {
                const amrex::Real term = std::sin(
                    std::numbers::pi_v<amrex::Real> * 0.5_rt *
                    (x - damping_west_start) /
                    (damping_west_end - damping_west_start));
                horizontal_coeff_west = term * term;
            }

            // North damping
            if (y < damping_north_start) {
                horizontal_coeff_north = 0.0_rt;
            } else if (y >= damping_north_end) {
                horizontal_coeff_north = 1.0_rt;
            } else {
                const amrex::Real term = std::sin(
                    std::numbers::pi_v<amrex::Real> * 0.5_rt *
                    (y - damping_north_start) /
                    (damping_north_end - damping_north_start));
                horizontal_coeff_north = term * term;
            }

            // South damping
            if (y > damping_south_start) {
                horizontal_coeff_south = 0.0_rt;
            } else if (y <= damping_south_end) {
                horizontal_coeff_south = 1.0_rt;
            } else {
                const amrex::Real term = std::sin(
                    std::numbers::pi_v<amrex::Real> * 0.5_rt *
                    (y - damping_south_start) /
                    (damping_south_end - damping_south_start));
                horizontal_coeff_south = term * term;
            }

            // Horizontal vertical coefficient
            if (z <= horizontal_abl_height) {
                vertical_coeff = 0.0_rt;
            } else if (z > z_sloped) {
                vertical_coeff = 1.0_rt;
            } else {
                const amrex::Real term = std::sin(
                    std::numbers::pi_v<amrex::Real> * 0.5_rt *
                    (z - horizontal_abl_height) /
                    (z_sloped - horizontal_abl_height));
                vertical_coeff = term * term;
            }

            levelDamping[nbx](i, j, k, 0) =
                vertical_coeff *
                (horizontal_coeff_east + horizontal_coeff_north +
                 horizontal_coeff_west + horizontal_coeff_south);

            // Add the full vertical damping
            if (z <= vertical_slope) {
                vertical_coeff = 0.0_rt;
            } else if (z > vertical_full) {
                vertical_coeff = 1.0_rt;
            } else {
                const amrex::Real term = std::sin(
                    std::numbers::pi_v<amrex::Real> * 0.5_rt *
                    (z - vertical_slope) /
                    (vertical_full - vertical_slope + 1e-15_rt));
                vertical_coeff = term * term;
            }

            levelDamping[nbx](i, j, k, 0) =
                amrex::min(
                    vertical_coeff + levelDamping[nbx](i, j, k, 0), 1.0_rt) /
                horizontal_tau;
        });
}

void PartialTerrainDrag::post_regrid_actions()
{
    const int nlevels = m_sim.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        initialize_fields(lev, m_sim.repo().mesh().Geom(lev));
    }
}

} // namespace kynema_sgf::partialterraindrag
