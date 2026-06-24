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
        } else if (m_type_forest == 2) {
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
            return {m_bbox_xlo, m_bbox_ylo, m_bbox_zlo,
                    m_bbox_xhi, m_bbox_yhi, m_bbox_zhi};
        }

        const amrex::Real search_tol = 1.1_rt;
        const amrex::Real search_radius =
            0.5_rt * search_tol * m_diameter_forest;
        const auto x0 = m_x_forest - search_radius;
        const auto y0 = m_y_forest - search_radius;
        const auto z0 = prob_lo[2];
        const auto x1 = m_x_forest + search_radius;
        const auto y1 = m_y_forest + search_radius;
        const auto z1 = m_height_forest;
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
};
} // namespace kynema_sgf::forestdrag

#endif
