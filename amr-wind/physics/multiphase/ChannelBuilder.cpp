#include "amr-wind/physics/multiphase/ChannelBuilder.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "amr-wind/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "amr-wind/utilities/IOManager.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::channelbuilder {

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
                // end Normal vector of plane: from start to end
                const amrex::Real a = end[0] - start[0];
                const amrex::Real b = end[1] - start[1];
                const amrex::Real c = end[2] - start[2];

                const amrex::Real p1 = a * (x - start[0]) + b * (y - start[1]) +
                                       c * (z - start[2]);
                const amrex::Real p2 =
                    a * (x - end[0]) + b * (y - end[1]) + c * (z - end[2]);
                // If point is within planes, check shape
                if (p1 * p2 <= 0.0_rt) {
                    // Convert coordinates to local segment coordinates
                    // Translate: x' = x - xstart, y' = y - ystart
                    // Rotate: x'' = x' * cos(theta) + y' * sin(theta)
                    //         y'' = x' * sin(theta) - y' * cos(theta)
                    // (this is a passive rotation; points are stationary and
                    // axes are rotated)
                    const amrex::Real xp = x - start[0];
                    const amrex::Real yp = y - start[1];
                    // Normalize normal vector to get xy rotation angle
                    const amrex::Real mag_xy = std::sqrt(a * a + b * b);
                    const amrex::Real cos_theta_xy = a / mag_xy;
                    const amrex::Real sin_theta_xy = b / mag_xy;
                    // For the moment, nothing is a function along the segment
                    // (fixed cross section)
                    const amrex::Real xpp =
                        xp * cos_theta_xy + yp * sin_theta_xy;
                    const amrex::Real ypp =
                        xp * sin_theta_xy - yp * cos_theta_xy;
                    // Now do the same with the vertical coordinate and xpp
                    const amrex::Real zp = z - start[2];
                    const amrex::Real mag = std::sqrt(a * a + b * b + c * c);
                    const amrex::Real cos_theta_xpz = mag_xy / mag;
                    const amrex::Real sin_theta_xpz = c / mag;
                    // Coordinates local to the segment centerline:
                    // x is along the segment, y is lateral/horizontal, and z is
                    // vertical to the segment
                    const amrex::Real xloc =
                        xpp * cos_theta_xpz + zp * sin_theta_xpz;
                    const amrex::Real zloc =
                        xpp * sin_theta_xpz - zp * cos_theta_xpz;
                    const amrex::Real yloc = ypp;

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
