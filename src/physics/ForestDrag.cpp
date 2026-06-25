#include "src/physics/ForestDrag.H"
#include "src/CFDSim.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParReduce.H"
#include "src/utilities/trig_ops.H"
#include "src/utilities/IOManager.H"
#include "AMReX_REAL.H"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <sstream>
#include <vector>
#include <ranges>

using namespace amrex::literals;

namespace kynema_sgf::forestdrag {

namespace {

// Strip inline comments and whitespace-only rows from point-cloud data files.
bool parse_data_line(const std::string& raw_line, std::istringstream& iss)
{
    const auto pos = raw_line.find('#');
    const auto line = raw_line.substr(0, pos);
    iss.clear();
    iss.str(line);
    return (
        !line.empty() &&
        line.find_first_not_of(" \t\r\n") != std::string::npos);
}

// Compute a 2D convex hull in the x-y plane using the monotonic chain method.
int compute_convex_hull_2d(
    const std::vector<std::pair<amrex::Real, amrex::Real>>& points,
    std::vector<std::pair<amrex::Real, amrex::Real>>& hull)
{
    if (points.size() < 3) {
        hull = points;
        return static_cast<int>(points.size());
    }

    std::vector<std::pair<amrex::Real, amrex::Real>> sorted_points = points;
    std::sort(
        sorted_points.begin(), sorted_points.end(),
        [](const auto& a, const auto& b) {
            return a.first < b.first ||
                   (std::abs(a.first - b.first) < constants::TIGHT_TOL &&
                    a.second < b.second);
        });

    auto cross = [](const auto& o, const auto& a, const auto& b) {
        return (a.first - o.first) * (b.second - o.second) -
               (a.second - o.second) * (b.first - o.first);
    };

    std::vector<std::pair<amrex::Real, amrex::Real>> lower;
    for (const auto& p : sorted_points) {
        while (lower.size() >= 2 &&
               cross(lower[lower.size() - 2], lower.back(), p) <=
                   constants::TIGHT_TOL) {
            lower.pop_back();
        }
        lower.push_back(p);
    }

    std::vector<std::pair<amrex::Real, amrex::Real>> upper;
    for (const auto& p : sorted_points | std::views::reverse) {
        while (upper.size() >= 2 &&
               cross(upper[upper.size() - 2], upper.back(), p) <=
                   constants::TIGHT_TOL) {
            upper.pop_back();
        }
        upper.push_back(p);
    }

    lower.pop_back();
    upper.pop_back();
    lower.insert(lower.end(), upper.begin(), upper.end());

    hull = lower;
    return static_cast<int>(hull.size());
}

} // namespace

ForestDrag::ForestDrag(CFDSim& sim)
    : m_sim(sim)
    , m_forest_drag(sim.repo().declare_field("forest_drag", 1, 1, 1))
    , m_forest_id(sim.repo().declare_field("forest_id", 1, 1, 1))
{

    // Accept either legacy cylinder input or point-cloud input, but not both.
    amrex::ParmParse pp(identifier());
    pp.query("forest_file", m_forest_file);
    const auto cyl_forest = pp.contains("forest_file");
    const auto point_forest = pp.contains("point_cloud_files");
    if (cyl_forest && point_forest) {
        amrex::Abort(
            "ForestDrag: Cannot specify both 'forest_file' and "
            "'point_cloud_files'");
    }
    pp.queryarr("point_cloud_files", m_point_cloud_files);
    if (point_forest) {
        // One drag coefficient is required per point-cloud forest file.
        pp.getarr("coefficients_of_drag", m_forest_cd);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            m_forest_cd.size() == m_point_cloud_files.size(),
            "ForestDrag: 'coefficients_of_drag' must have the same "
            "number of entries as 'point_cloud_files'");
        pp.query("point_neighbors", m_point_neighbors);
        pp.query("point_interp_eps", m_point_interp_eps);
        // Keep neighbor count bounded by fixed-size device arrays in kernels.
        m_point_neighbors =
            amrex::max<int>(1, amrex::min<int>(m_point_neighbors, 8));
    }

    // Register outputs and initialize field defaults.
    m_sim.io_manager().register_output_var("forest_drag");
    m_sim.io_manager().register_output_var("forest_id");

    m_forest_drag.setVal(0.0_rt);
    m_forest_id.setVal(-1.0_rt);
    m_forest_id.set_default_fillpatch_bc(m_sim.time());
    m_forest_drag.set_default_fillpatch_bc(m_sim.time());
}

void ForestDrag::initialize_fields(int level, const amrex::Geometry& geom)
{
    BL_PROFILE("kynema-sgf::" + this->identifier() + "::initialize_fields");

    // Build host-side forest metadata for the requested AMR level.
    amrex::Vector<ForestPoint> cloud_points;
    amrex::Vector<kynema_sgf::forestdrag::ForestHullEdge> hull_edges;
    amrex::Vector<Forest> forests;
    if (!m_point_cloud_files.empty()) {
        forests = read_point_cloud_forests(level, cloud_points, hull_edges);
    } else {
        forests = read_cylinder_forests(level);
    }

    // Mirror host vectors on device so kernels can access forest descriptors,
    // point samples, and precomputed hull half-space edges.
    amrex::Gpu::DeviceVector<Forest> d_forests(forests.size());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, forests.begin(), forests.end(),
        d_forests.begin());
    amrex::Gpu::DeviceVector<ForestPoint> d_cloud_points(cloud_points.size());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, cloud_points.begin(), cloud_points.end(),
        d_cloud_points.begin());

    amrex::Gpu::DeviceVector<kynema_sgf::forestdrag::ForestHullEdge>
        d_hull_edges(hull_edges.size());
    if (!hull_edges.empty()) {
        amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, hull_edges.begin(), hull_edges.end(),
            d_hull_edges.begin());
    }

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    auto& drag = m_forest_drag(level);
    auto& fst_id = m_forest_id(level);

    // Rebuild fields from scratch each time we initialize this level.
    drag.setVal(0.0_rt);
    fst_id.setVal(-1.0_rt);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(m_forest_drag(level)); mfi.isValid(); ++mfi) {
        const auto& vbx = mfi.growntilebox();
        for (int nf = 0; nf < static_cast<int>(forests.size()); nf++) {
            // First trim work to a forest-local bounding box before per-cell
            // logic.
            const auto bxi = vbx & forests[nf].bounding_box(geom);
            if (!bxi.isEmpty()) {
                const auto& levelDrag = drag.array(mfi);
                const auto& levelId = fst_id.array(mfi);
                const auto* forests_ptr = d_forests.data();
                const auto* points_ptr = d_cloud_points.data();
                const auto* hull_edges_ptr = d_hull_edges.data();
                const int num_neighbors = m_point_neighbors;
                const amrex::Real interp_eps = m_point_interp_eps;
                amrex::ParallelFor(
                    bxi, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        // Convert integer cell indices to cell-center
                        // coordinates.
                        const auto x = prob_lo[0] + ((i + 0.5_rt) * dx[0]);
                        const auto y = prob_lo[1] + ((j + 0.5_rt) * dx[1]);
                        const auto z = prob_lo[2] + ((k + 0.5_rt) * dx[2]);
                        const auto& fst = forests_ptr[nf];

                        // Legacy mode: cylindrical footprint with vertical LAD
                        // profile.
                        if (fst.m_drag_mode == 0) {
                            const auto radius = std::sqrt(
                                ((x - fst.m_x_forest) * (x - fst.m_x_forest)) +
                                ((y - fst.m_y_forest) * (y - fst.m_y_forest)));
                            if (z <= fst.m_height_forest &&
                                radius <= (0.5_rt * fst.m_diameter_forest)) {
                                const auto treelaimax = fst.lm();
                                levelId(i, j, k) = fst.m_id;
                                levelDrag(i, j, k) +=
                                    fst.m_cd_forest *
                                    fst.area_fraction(z, treelaimax);
                            }
                        } else if (
                            fst.m_cloud_point_count > 0 &&
                            fst.point_in_on_hull(x, y, hull_edges_ptr)) {
                            // Point-cloud mode: only evaluate interpolation for
                            // cells inside the convex hull projected in x-y.
                            constexpr int max_neighbors = 8;
                            constexpr amrex::Real huge = constants::LARGE_NUM;
                            amrex::Array<amrex::Real, max_neighbors> nearest_d2;
                            amrex::Array<amrex::Real, max_neighbors>
                                nearest_lad;
                            amrex::Array<amrex::Real, max_neighbors> nearest_z;

                            for (int n = 0; n < max_neighbors; ++n) {
                                nearest_d2[n] = huge;
                                nearest_lad[n] = 0.0_rt;
                                nearest_z[n] = 0.0_rt;
                            }

                            const int pstart = fst.m_cloud_point_offset;
                            const int pend = pstart + fst.m_cloud_point_count;
                            for (int p = pstart; p < pend; ++p) {
                                const auto& pt = points_ptr[p];
                                const auto dxp = x - pt.m_x;
                                const auto dyp = y - pt.m_y;
                                const auto dzp = z - pt.m_z;
                                const auto d2 =
                                    (dxp * dxp) + (dyp * dyp) + (dzp * dzp);

                                // Maintain a sorted nearest-neighbor list
                                // in-place.
                                if (d2 < nearest_d2[num_neighbors - 1]) {
                                    int insert = num_neighbors - 1;
                                    while (insert > 0 &&
                                           d2 < nearest_d2[insert - 1]) {
                                        nearest_d2[insert] =
                                            nearest_d2[insert - 1];
                                        nearest_lad[insert] =
                                            nearest_lad[insert - 1];
                                        nearest_z[insert] =
                                            nearest_z[insert - 1];
                                        --insert;
                                    }
                                    nearest_d2[insert] = d2;
                                    nearest_lad[insert] = pt.m_lad;
                                    nearest_z[insert] = pt.m_z;
                                }
                            }

                            // Limit interpolation above canopy by nearest
                            // sample heights.
                            amrex::Real max_z_neighbors = 0.0_rt;
                            for (int n = 0; n < num_neighbors; ++n) {
                                max_z_neighbors = amrex::max<amrex::Real>(
                                    max_z_neighbors, nearest_z[n]);
                            }

                            // Interpolate LAD with inverse-distance weighting.
                            // If we are essentially on a sample point, use it
                            // directly.
                            const auto eps2 = interp_eps * interp_eps;
                            amrex::Real lad_interp = 0.0_rt;
                            if (nearest_d2[0] < eps2) {
                                lad_interp = nearest_lad[0];
                            } else if (z - 0.5_rt * dx[2] <= max_z_neighbors) {
                                amrex::Real sum_w = 0.0_rt;
                                amrex::Real sum_lad = 0.0_rt;
                                for (int n = 0; n < num_neighbors; ++n) {
                                    if (nearest_d2[n] >= huge) {
                                        continue;
                                    }
                                    const auto w =
                                        1.0_rt /
                                        std::sqrt(nearest_d2[n] + eps2);
                                    sum_w += w;
                                    sum_lad += w * nearest_lad[n];
                                }
                                if (sum_w > 0.0_rt) {
                                    lad_interp = sum_lad / sum_w;
                                }
                            }

                            // Apply drag contribution only for positive
                            // interpolated LAD.
                            if (lad_interp > 0.0_rt) {
                                levelId(i, j, k) = fst.m_id;
                                levelDrag(i, j, k) +=
                                    fst.m_cd_forest * lad_interp;
                            }
                        }
                    });
            }
        }
    }
}

void ForestDrag::post_regrid_actions()
{
    // Forest fields depend on geometry and valid boxes, so recompute on regrid.
    const int nlevels = m_sim.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        initialize_fields(lev, m_sim.repo().mesh().Geom(lev));
    }
}

amrex::Vector<Forest> ForestDrag::read_cylinder_forests(const int level) const
{
    BL_PROFILE("kynema-sgf::" + this->identifier() + "::read_cylinder_forests");

    // Legacy format: one row per cylindrical forest.
    std::ifstream file(m_forest_file, std::ios::in);
    if (!file.good()) {
        amrex::Abort("Cannot find file " + m_forest_file);
    }

    //! TreeType xc yc height diameter cd lai laimax
    amrex::Vector<Forest> forests;
    const auto& geom = m_sim.repo().mesh().Geom(level);
    const auto& ba = m_sim.repo().mesh().boxArray(level);
    int cnt = 0;
    amrex::Real value1, value2, value3, value4, value5, value6, value7, value8;
    while (file >> value1 >> value2 >> value3 >> value4 >> value5 >> value6 >>
           value7 >> value8) {
        Forest f;
        f.m_id = cnt;
        f.m_type_forest = value1;
        f.m_x_forest = value2;
        f.m_y_forest = value3;
        f.m_height_forest = value4;
        f.m_diameter_forest = value5;
        f.m_cd_forest = value6;
        f.m_lai_forest = value7;
        f.m_laimax_forest = value8;

        // Keep only forests intersecting this rank's local level boxes.
        const auto bx = f.bounding_box(geom);
        if (ba.intersects(bx)) {
            forests.push_back(f);
        }
        cnt++;
    }
    file.close();
    return forests;
}

amrex::Vector<Forest> ForestDrag::read_point_cloud_forests(
    const int level,
    amrex::Vector<ForestPoint>& points,
    amrex::Vector<kynema_sgf::forestdrag::ForestHullEdge>& hull_edges) const
{
    BL_PROFILE(
        "kynema-sgf::" + this->identifier() + "::read_point_cloud_forests");

    // Point-cloud mode builds one Forest descriptor per input file.
    amrex::Vector<Forest> forests;
    const auto& geom = m_sim.repo().mesh().Geom(level);
    const auto& ba = m_sim.repo().mesh().boxArray(level);

    int cnt = 0;
    for (int i = 0; i < m_point_cloud_files.size(); ++i) {

        const auto& cloud_file = m_point_cloud_files[i];

        std::ifstream cloud_data(cloud_file, std::ios::in);
        if (!cloud_data.good()) {
            amrex::Abort("Cannot find forest point cloud file " + cloud_file);
        }

        Forest f;
        f.m_id = cnt;
        f.m_drag_mode = 1;
        f.m_cd_forest = m_forest_cd[i];
        f.m_cloud_point_offset = static_cast<int>(points.size());

        // Track extents for a quick coarse bounding box test.
        amrex::Real xmin = std::numeric_limits<amrex::Real>::max();
        amrex::Real ymin = std::numeric_limits<amrex::Real>::max();
        amrex::Real zmin = std::numeric_limits<amrex::Real>::max();
        amrex::Real xmax = std::numeric_limits<amrex::Real>::lowest();
        amrex::Real ymax = std::numeric_limits<amrex::Real>::lowest();
        amrex::Real zmax = std::numeric_limits<amrex::Real>::lowest();

        std::vector<std::pair<amrex::Real, amrex::Real>> xy_points;

        std::string p_line;
        std::istringstream p_stream;
        while (std::getline(cloud_data, p_line)) {
            if (!parse_data_line(p_line, p_stream)) {
                continue;
            }

            ForestPoint pt;
            if (!(p_stream >> pt.m_x >> pt.m_y >> pt.m_z >> pt.m_lad)) {
                amrex::Abort(
                    "Malformed point entry in " + cloud_file +
                    ". Expected: x y z lad");
            }

            points.push_back(pt);
            xy_points.emplace_back(pt.m_x, pt.m_y);
            xmin = amrex::min<amrex::Real>(xmin, pt.m_x);
            ymin = amrex::min<amrex::Real>(ymin, pt.m_y);
            zmin = amrex::min<amrex::Real>(zmin, pt.m_z);
            xmax = amrex::max<amrex::Real>(xmax, pt.m_x);
            ymax = amrex::max<amrex::Real>(ymax, pt.m_y);
            zmax = amrex::max<amrex::Real>(zmax, pt.m_z);
        }

        // Number of points belonging to this forest in the flattened point
        // array.
        f.m_cloud_point_count =
            static_cast<int>(points.size()) - f.m_cloud_point_offset;
        if (f.m_cloud_point_count <= 0) {
            amrex::Abort(
                "Forest point cloud file has no valid points: " + cloud_file);
        }

        // Precompute the x-y convex hull used as an interpolation-domain
        // filter, then convert it to inward half-space edges.
        std::vector<std::pair<amrex::Real, amrex::Real>> hull_2d;
        compute_convex_hull_2d(xy_points, hull_2d);
        if (hull_2d.size() < 3) {
            amrex::Abort(
                "Forest point cloud file must contain at least 3 non-collinear "
                "x-y points: " +
                cloud_file);
        }

        // Ensure CCW winding so inward normals are consistently defined.
        amrex::Real twice_area = 0.0_rt;
        for (int v = 0; v < static_cast<int>(hull_2d.size()); ++v) {
            const int vn = (v + 1) % static_cast<int>(hull_2d.size());
            twice_area += (hull_2d[v].first * hull_2d[vn].second) -
                          (hull_2d[vn].first * hull_2d[v].second);
        }
        if (twice_area < 0.0_rt) {
            std::reverse(hull_2d.begin(), hull_2d.end());
        }

        f.m_hull_edge_offset = static_cast<int>(hull_edges.size());
        f.m_hull_edge_count = static_cast<int>(hull_2d.size());
        for (int v = 0; v < static_cast<int>(hull_2d.size()); ++v) {
            const int vn = (v + 1) % static_cast<int>(hull_2d.size());
            const auto xi = hull_2d[v].first;
            const auto yi = hull_2d[v].second;
            const auto xj = hull_2d[vn].first;
            const auto yj = hull_2d[vn].second;
            const auto ex = xj - xi;
            const auto ey = yj - yi;
            const auto nx = -ey;
            const auto ny = ex;
            hull_edges.emplace_back(
                ForestHullEdge{
                    nx, ny, (nx * xi) + (ny * yi),
                    amrex::Math::abs(nx) + amrex::Math::abs(ny) + 1.0_rt});
        }

        // Slightly pad the point-cloud extents to avoid missing edge cells due
        // to floating-point and cell-center alignment effects.
        const amrex::Real pad =
            0.5_rt *
            (geom.CellSizeArray()[0] + geom.CellSizeArray()[1] +
             geom.CellSizeArray()[2]) /
            3.0_rt;
        f.m_bbox_xlo = xmin - pad;
        f.m_bbox_ylo = ymin - pad;
        f.m_bbox_zlo = zmin - pad;
        f.m_bbox_xhi = xmax + pad;
        f.m_bbox_yhi = ymax + pad;
        f.m_bbox_zhi = zmax + pad;

        // Keep legacy fields coherent for diagnostics/output.
        f.m_x_forest = 0.5_rt * (xmin + xmax);
        f.m_y_forest = 0.5_rt * (ymin + ymax);
        f.m_height_forest = zmax;
        f.m_diameter_forest =
            2.0_rt *
            amrex::max<amrex::Real>(xmax - f.m_x_forest, ymax - f.m_y_forest);
        f.m_type_forest = 0.0_rt;
        f.m_lai_forest = 0.0_rt;
        f.m_laimax_forest = 0.0_rt;

        // Keep only forests intersecting this rank's local level boxes.
        if (ba.intersects(f.bounding_box(geom))) {
            forests.push_back(f);
        }
        ++cnt;
    }

    return forests;
}
} // namespace kynema_sgf::forestdrag
