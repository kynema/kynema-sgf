#ifndef FORESTDRAG_H
#define FORESTDRAG_H

#include "src/core/Physics.H"
#include "src/core/Field.H"
#include "src/CFDSim.H"
#include "src/utilities/index_operations.H"
#include "src/utilities/integrals.H"
#include "src/utilities/constants.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::forestdrag {

// One precomputed inward-facing half-space edge for convex hull tests.
struct ForestHullEdge
{
    amrex::Real m_nx;
    amrex::Real m_ny;
    amrex::Real m_d;
    amrex::Real m_tol_scale;
};

// Per-forest descriptor consumed by host logic and device kernels.
struct Forest
{
    // Forest identifier written to forest_id field.
    int m_id{-1};
    // 0 = legacy cylindrical forest, 1 = point-cloud forest.
    int m_drag_mode{0};

    // Point-cloud storage slice [offset, offset + count) in flattened point
    // array.
    int m_cloud_point_offset{0};
    int m_cloud_point_count{0};

    // Convex hull half-space storage slice in flattened edge array.
    int m_hull_edge_offset{0};
    int m_hull_edge_count{0};

    // Legacy cylindrical/model parameters.
    amrex::Real m_type_forest;
    amrex::Real m_x_forest;
    amrex::Real m_y_forest;
    amrex::Real m_height_forest;
    amrex::Real m_diameter_forest;
    amrex::Real m_cd_forest;
    amrex::Real m_lai_forest;
    amrex::Real m_laimax_forest;

    // Axis-aligned bounding box used for coarse spatial pruning.
    amrex::Real m_bbox_xlo{0.0_rt};
    amrex::Real m_bbox_ylo{0.0_rt};
    amrex::Real m_bbox_zlo{0.0_rt};
    amrex::Real m_bbox_xhi{0.0_rt};
    amrex::Real m_bbox_yhi{0.0_rt};
    amrex::Real m_bbox_zhi{0.0_rt};

    // Range of the terrain height on this level. Both stay zero without
    // TerrainDrag, so the bounding boxes are unchanged in that case.
    amrex::Real m_terrain_zmin{0.0_rt};
    amrex::Real m_terrain_zmax{0.0_rt};

    // Compute the type-2 peak LAD scaling factor from integrated LAI.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real lm() const
    {
        amrex::Real treelaimax = 0.0_rt;
        if (m_type_forest == 2) {
            const amrex::Real captured_height_forest = m_height_forest;
            const amrex::Real treeZm = m_laimax_forest * captured_height_forest;
            const int n = 100;
            const amrex::Real int_0_fh = kynema_sgf::utils::trapz(
                0, captured_height_forest - constants::TIGHT_TOL, n,
                [=](const amrex::Real x) {
                    const amrex::Real ratio =
                        (captured_height_forest - treeZm) /
                        (captured_height_forest - x);
                    const auto exponent = (x < treeZm) ? 6.0_rt : 0.5_rt;
                    return std::pow(ratio, exponent) *
                           std::exp(exponent * (1 - ratio));
                });
            treelaimax = m_lai_forest / int_0_fh;
        }
        return treelaimax;
    }

    // Vertical LAD profile used by legacy cylindrical forests.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    area_fraction(const amrex::Real z, const amrex::Real treelaimax) const
    {
        amrex::Real af = 0.0_rt;
        if (m_type_forest == 1) {
            af = m_lai_forest / m_height_forest;
        } else if (m_type_forest == 2 && z < m_height_forest) {
            // The profile tends to zero at the tree top, where the ratio
            // below divides by zero, so the top and above are left at zero.
            const auto treeZm = m_laimax_forest * m_height_forest;
            const auto ratio =
                (m_height_forest - treeZm) / (m_height_forest - z);
            const auto exponent = (z < treeZm) ? 6.0_rt : 0.5_rt;
            af = treelaimax * std::pow(ratio, exponent) *
                 std::exp(exponent * (1 - ratio));
        }
        return af;
    }

    // Convex-hull half-space test (boundary-inclusive) using precomputed
    // inward edge normals.
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool point_in_on_hull(
        const amrex::Real px,
        const amrex::Real py,
        const ForestHullEdge* hull_edges) const
    {
        if (m_hull_edge_count <= 0) {
            return true;
        }
        for (int i = 0; i < m_hull_edge_count; ++i) {
            const auto& e = hull_edges[m_hull_edge_offset + i];
            const amrex::Real signed_dist =
                (e.m_nx * px) + (e.m_ny * py) - e.m_d;
            if (signed_dist < -constants::TIGHT_TOL * e.m_tol_scale) {
                return false;
            }
        }
        return true;
    }

    // Return real-space bounds for either point-cloud or cylindrical mode.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::RealBox real_bounding_box(
        const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& prob_lo) const
    {
        if (m_drag_mode == 1) {
            return {m_bbox_xlo, m_bbox_ylo, m_bbox_zlo + m_terrain_zmin,
                    m_bbox_xhi, m_bbox_yhi, m_bbox_zhi + m_terrain_zmax};
        }

        const amrex::Real search_tol = 1.1_rt;
        const amrex::Real search_radius =
            0.5_rt * search_tol * m_diameter_forest;
        const auto x0 = m_x_forest - search_radius;
        const auto y0 = m_y_forest - search_radius;
        const auto z0 = prob_lo[2];
        const auto x1 = m_x_forest + search_radius;
        const auto y1 = m_y_forest + search_radius;
        const auto z1 = m_height_forest + m_terrain_zmax;
        return {x0, y0, z0, x1, y1, z1};
    }

    // Integer box used for tile intersection tests before launching kernels.
    [[nodiscard]] amrex::Box bounding_box(const amrex::Geometry& geom) const
    {
        return utils::realbox_to_box(
            real_bounding_box(geom.ProbLoArray()), geom);
    }
};

// A single LAD sample for point-cloud forests.
struct ForestPoint
{
    amrex::Real m_x;
    amrex::Real m_y;
    amrex::Real m_z;
    amrex::Real m_lad;
};

/** LAD of a legacy cylindrical forest at a point
 *
 *  \param fst Forest descriptor (m_drag_mode == 0)
 *  \param x Point x coordinate
 *  \param y Point y coordinate
 *  \param zagl Point height above the local ground
 *  \param lad Leaf area density at the point (output)
 *  \return true when the point lies inside the forest cylinder
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE bool cylinder_lad(
    const Forest& fst,
    const amrex::Real x,
    const amrex::Real y,
    const amrex::Real zagl,
    amrex::Real& lad)
{
    lad = 0.0_rt;
    const auto radius = std::sqrt(
        ((x - fst.m_x_forest) * (x - fst.m_x_forest)) +
        ((y - fst.m_y_forest) * (y - fst.m_y_forest)));
    if (zagl <= fst.m_height_forest &&
        radius <= (0.5_rt * fst.m_diameter_forest)) {
        lad = fst.area_fraction(zagl, fst.lm());
        return true;
    }
    return false;
}

/** LAD of a point-cloud forest at a point by inverse-distance weighting
 *
 *  \param fst Forest descriptor (m_drag_mode == 1)
 *  \param points Flattened point-cloud samples of all forests
 *  \param hull_edges Flattened convex hull edges of all forests
 *  \param x Point x coordinate
 *  \param y Point y coordinate
 *  \param z Point height, in the same frame as the sample heights
 *  \param half_dz Half the cell height, for the canopy-top cutoff
 *  \param num_neighbors Number of nearest samples to weight (at most 8)
 *  \param interp_eps Regularization length of the weights
 *  \return Interpolated LAD, zero outside the hull or above the canopy
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real point_cloud_lad(
    const Forest& fst,
    const ForestPoint* points,
    const ForestHullEdge* hull_edges,
    const amrex::Real x,
    const amrex::Real y,
    const amrex::Real z,
    const amrex::Real half_dz,
    const int num_neighbors,
    const amrex::Real interp_eps)
{
    if (fst.m_cloud_point_count <= 0 ||
        !fst.point_in_on_hull(x, y, hull_edges)) {
        return 0.0_rt;
    }

    const int max_neighbors = 8;
    const amrex::Real huge = constants::LARGE_NUM;
    amrex::GpuArray<amrex::Real, 8> nearest_d2;
    amrex::GpuArray<amrex::Real, 8> nearest_lad;
    amrex::GpuArray<amrex::Real, 8> nearest_z;
    for (int n = 0; n < max_neighbors; ++n) {
        nearest_d2[n] = huge;
        nearest_lad[n] = 0.0_rt;
        nearest_z[n] = 0.0_rt;
    }

    // Sorted nearest-neighbor list, updated in place.
    const int pstart = fst.m_cloud_point_offset;
    const int pend = pstart + fst.m_cloud_point_count;
    for (int p = pstart; p < pend; ++p) {
        const auto& pt = points[p];
        const auto dxp = x - pt.m_x;
        const auto dyp = y - pt.m_y;
        const auto dzp = z - pt.m_z;
        const auto d2 = (dxp * dxp) + (dyp * dyp) + (dzp * dzp);
        if (d2 < nearest_d2[num_neighbors - 1]) {
            int insert = num_neighbors - 1;
            while (insert > 0 && d2 < nearest_d2[insert - 1]) {
                nearest_d2[insert] = nearest_d2[insert - 1];
                nearest_lad[insert] = nearest_lad[insert - 1];
                nearest_z[insert] = nearest_z[insert - 1];
                --insert;
            }
            nearest_d2[insert] = d2;
            nearest_lad[insert] = pt.m_lad;
            nearest_z[insert] = pt.m_z;
        }
    }

    // No LAD above the highest of the nearest samples.
    amrex::Real max_z_neighbors = 0.0_rt;
    for (int n = 0; n < num_neighbors; ++n) {
        max_z_neighbors =
            amrex::max<amrex::Real>(max_z_neighbors, nearest_z[n]);
    }

    const auto eps2 = interp_eps * interp_eps;
    amrex::Real lad_interp = 0.0_rt;
    if (nearest_d2[0] < eps2) {
        lad_interp = nearest_lad[0];
    } else if (z - half_dz <= max_z_neighbors) {
        amrex::Real sum_w = 0.0_rt;
        amrex::Real sum_lad = 0.0_rt;
        for (int n = 0; n < num_neighbors; ++n) {
            if (nearest_d2[n] >= huge) {
                continue;
            }
            const auto w = 1.0_rt / std::sqrt(nearest_d2[n] + eps2);
            sum_w += w;
            sum_lad += w * nearest_lad[n];
        }
        if (sum_w > 0.0_rt) {
            lad_interp = sum_lad / sum_w;
        }
    }
    return lad_interp;
}

/** Forestdrag Flow physics
 *  \ingroup physics
 */

class ForestDrag : public Physics::Register<ForestDrag>
{
public:
    static std::string identifier() { return "ForestDrag"; }

    explicit ForestDrag(CFDSim& sim);

    ~ForestDrag() override = default;

    void
    initialize_fields(int /*level*/, const amrex::Geometry& /*geom*/) override;

    void pre_init_actions() override {}

    void post_init_actions() override {}

    void post_regrid_actions() override;

    void pre_advance_work() override {}

    void post_advance_work() override {}

    // Parse legacy cylindrical forest file for a specific AMR level.
    [[nodiscard]] amrex::Vector<Forest> read_cylinder_forests(int level) const;

    // Parse point-cloud forests, flatten points/hull-edges, and build forest
    // descriptors.
    [[nodiscard]] amrex::Vector<Forest> read_point_cloud_forests(
        int level,
        amrex::Vector<ForestPoint>& points,
        amrex::Vector<ForestHullEdge>& hull_edges) const;

private:
    CFDSim& m_sim;
    // Output fields populated by this physics source.
    Field& m_forest_drag;
    Field& m_forest_id;

    // Input configuration for legacy and point-cloud modes.
    std::string m_forest_file{"forest.amrwind"};
    amrex::Vector<std::string> m_point_cloud_files;
    amrex::Vector<amrex::Real> m_forest_cd;

    // Point-cloud interpolation controls.
    int m_point_neighbors{4};
    amrex::Real m_point_interp_eps{constants::TIGHT_TOL};

    // Place the forests on the TerrainDrag terrain: heights in the forest and
    // point-cloud files are taken above the local ground.
    bool m_terrain_aware{true};
    // TerrainDrag was created before ForestDrag, so its terrain_height field
    // is filled first on every (re)initialization.
    bool m_terrain_created_first{false};
    // Terrain height range on the level being initialized, used to extend the
    // forest bounding boxes.
    amrex::Real m_terrain_zmin{0.0_rt};
    amrex::Real m_terrain_zmax{0.0_rt};
};
} // namespace kynema_sgf::forestdrag

#endif
