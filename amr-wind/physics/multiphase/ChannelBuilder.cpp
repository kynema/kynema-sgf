#include "amr-wind/physics/multiphase/ChannelBuilder.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/IOManager.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::channelbuilder {

[[nodiscard]] bool trapezoid(
    const amrex::Real& top,
    const amrex::Real& bottom,
    const amrex::Real& height,
    const amrex::Real& hcoord,
    const amrex::Real& vcoord)
{
    return (vcoord >= -height / 2.0_rt) && (vcoord <= height / 2.0_rt) &&
           (hcoord >=
            (-(top + bottom) / 2.0_rt + ((bottom - top) / height) * vcoord)) &&
           (hcoord <=
            ((top + bottom) / 2.0_rt - ((bottom - top) / height) * vcoord));
}

[[nodiscard]] bool ellipse(
    const amrex::Real& ax_horz,
    const amrex::Real& ax_vert,
    const amrex::Real& hcoord,
    const amrex::Real& vcoord)
{
    return (
        (hcoord * hcoord) / (ax_horz * ax_horz) +
            (vcoord * vcoord) / (ax_vert * ax_vert) <=
        1.0_rt);
}

[[nodiscard]] AMREX_GPU_HOST_DEVICE bool is_point_within_planes(
    const amrex::Real& x,
    const amrex::Real& y,
    const amrex::Real& z,
    const amrex::Real& start_x,
    const amrex::Real& start_y,
    const amrex::Real& start_z,
    const amrex::Real& end_x,
    const amrex::Real& end_y,
    const amrex::Real& end_z)
{
    // Normal vector of plane: from start to end
    const amrex::Real a = end_x - start_x;
    const amrex::Real b = end_y - start_y;
    const amrex::Real c = end_z - start_z;

    // Dot product with plane normal at start and end
    const amrex::Real p1 =
        a * (x - start_x) + b * (y - start_y) + c * (z - start_z);
    const amrex::Real p2 = a * (x - end_x) + b * (y - end_y) + c * (z - end_z);

    // Point is within planes if it's on opposite sides or on a plane
    return (p1 * p2 <= 0.0_rt);
}

[[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::GpuArray<amrex::Real, 3>
transform_to_local_coordinates(
    const amrex::Real& x,
    const amrex::Real& y,
    const amrex::Real& z,
    const amrex::Real& start_x,
    const amrex::Real& start_y,
    const amrex::Real& start_z,
    const amrex::Real& end_x,
    const amrex::Real& end_y,
    const amrex::Real& end_z)
{
    // Segment direction vector
    const amrex::Real a = end_x - start_x;
    const amrex::Real b = end_y - start_y;
    const amrex::Real c = end_z - start_z;

    // Translate point relative to segment start
    const amrex::Real xp = x - start_x;
    const amrex::Real yp = y - start_y;
    const amrex::Real zp = z - start_z;

    // Rotate around z-axis based on xy component of direction
    const amrex::Real mag_xy = std::sqrt(a * a + b * b);
    const amrex::Real cos_theta_xy = a / mag_xy;
    const amrex::Real sin_theta_xy = b / mag_xy;

    const amrex::Real xpp =
        xp * cos_theta_xy + yp * sin_theta_xy;
    const amrex::Real ypp =
        - xp * sin_theta_xy + yp * cos_theta_xy;

    // Rotate around y-axis based on z component
    const amrex::Real mag = std::sqrt(a * a + b * b + c * c);
    const amrex::Real cos_theta_xpz = mag_xy / mag;
    const amrex::Real sin_theta_xpz = c / mag;

    // Local coordinates: xloc along segment, yloc lateral, zloc vertical
    const amrex::Real xloc =
        xpp * cos_theta_xpz + zp * sin_theta_xpz;
    const amrex::Real zloc =
        - xpp * sin_theta_xpz + zp * cos_theta_xpz;

    return amrex::GpuArray<amrex::Real, 3>{
        {xloc, ypp, zloc}};
}

ChannelBuilder::ChannelBuilder(CFDSim& sim)
    : m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_terrain_blank(sim.repo().declare_int_field("terrain_blank", 1, 1, 1))
{

    m_sim.io_manager().register_output_int_var("terrain_blank");

    m_terrain_blank.setVal(0);
    amrex::Vector<std::string> labels;
    amrex::ParmParse pp(identifier());
    is_multiphase = pp.contains("water_level");
    if (is_multiphase) {
        pp.get("water_level", m_water_level);
        pp.get("land_level", m_land_level);
        amrex::ParmParse pp_multiphase("MultiPhase");
        pp_multiphase.add("water_level", m_water_level);
    }
    pp.getarr("segment_labels", labels);
    for (const auto& lbl : labels) {
        const std::string key = identifier() + "." + lbl;
        amrex::ParmParse pp1(key);
        std::string stype = "Ellipse";
        ChannelSegmentType type = ChannelSegmentType::Ellipse;
        amrex::Real dim0 = 0.0_rt;
        amrex::Real dim1 = 0.0_rt;
        amrex::Real dim2 = 0.0_rt;
        amrex::Vector<amrex::Real> seg_start{0.0_rt, 0.0_rt, 0.0_rt};
        amrex::Vector<amrex::Real> seg_end{0.0_rt, 0.0_rt, 0.0_rt};

        pp1.query("type", stype);
        if (stype == "Ellipse") {
            type = ChannelSegmentType::Ellipse;
            if (pp1.contains("diameter")) {
                pp1.get("diameter", dim0);
                dim1 = dim0;
            } else {
                pp1.get("horizontal_axis", dim0);
                pp1.get("vertical_axis", dim1);
            }
        } else if (stype == "Trapezoid") {
            type = ChannelSegmentType::Trapezoid;
            pp1.get("top_width", dim1);
            pp1.get("bottom_width", dim0);
            pp1.get("height", dim2);
        } else {
            amrex::Abort(
                "Invalid channel segment type specified: " + stype +
                ". Only 'Ellipse' and 'Trapezoid' are supported.");
        }
        pp1.getarr("segment_start", seg_start);
        pp1.getarr("segment_end", seg_end);

        m_type.emplace_back(type);
        m_dim0.emplace_back(dim0);
        m_dim1.emplace_back(dim1);
        m_dim2.emplace_back(dim2);
        m_segment_start.emplace_back(
            amrex::Array<amrex::Real, 3>{
                seg_start[0], seg_start[1], seg_start[2]});
        m_segment_end.emplace_back(
            amrex::Array<amrex::Real, 3>{seg_end[0], seg_end[1], seg_end[2]});
    }
}

void ChannelBuilder::initialize_fields(int level, const amrex::Geometry& geom)
{
    // if multiphase, check for multiphase physics
    if (is_multiphase && !m_sim.physics_manager().contains("MultiPhase")) {
        amrex::Abort(
            "ChannelBuilder: MultiPhase physics must be enabled to use "
            "multiphase channel builder");
    }

    const auto& dx = geom.CellSizeArray();
    const auto& prob_lo = geom.ProbLoArray();
    const auto& prob_hi = geom.ProbHiArray();
    auto& blank_mfab = m_terrain_blank(level);
    auto blank_arrs = blank_mfab.arrays();
    amrex::ParallelFor(
        blank_mfab, m_terrain_blank.num_grow(),
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) {
            const amrex::Real x = prob_lo[0] + ((i + 0.5_rt) * dx[0]);
            const amrex::Real y = prob_lo[1] + ((j + 0.5_rt) * dx[1]);
            const amrex::Real z = prob_lo[2] + ((k + 0.5_rt) * dx[2]);

            // Start with cell blanked
            bool outside_channel = true;

            // Loop through segments and determine if cell is within any channel
            // segment
            for (int seg = 0; seg < m_type.size(); ++seg) {
                const auto& start = m_segment_start[seg];
                const auto& end = m_segment_end[seg];
                const auto& seg_type = m_type[seg];
                const auto& dim0 = m_dim0[seg];
                const auto& dim1 = m_dim1[seg];
                const auto& dim2 = m_dim2[seg];

                // Check if point is within bounding planes of segment start and
                // end
                if (is_point_within_planes(
                        x, y, z, start[0], start[1], start[2], end[0], end[1],
                        end[2])) {
                    // Transform to local segment coordinates
                    const auto local_coords = transform_to_local_coordinates(
                        x, y, z, start[0], start[1], start[2], end[0], end[1],
                        end[2]);
                    const amrex::Real yloc = local_coords[1];
                    const amrex::Real zloc = local_coords[2];

                    if (seg_type == ChannelSegmentType::Ellipse) {
                        outside_channel &= !ellipse(dim0, dim1, yloc, zloc);
                    } else if (seg_type == ChannelSegmentType::Trapezoid) {
                        outside_channel &=
                            !trapezoid(dim0, dim1, dim2, yloc, zloc);
                    }
                }
            }

            // Do adjustments for multiphase case
            if (is_multiphase && z > m_land_level) {
                // Above land level means unblanked
                outside_channel = false;
            }

            blank_arrs[nbx](i, j, k) = static_cast<int>(outside_channel);
        });
    amrex::Gpu::streamSynchronize();

    // Do not set "drag" cells until improving drag forcing to handle different
    // directions (i.e., not just above terrain)

    // Same goes for roughness
}

void ChannelBuilder::post_regrid_actions()
{
    const int nlevels = m_sim.repo().num_active_levels();
    for (int lev = 0; lev < nlevels; ++lev) {
        initialize_fields(lev, m_sim.repo().mesh().Geom(lev));
    }
}

} // namespace amr_wind::channelbuilder
