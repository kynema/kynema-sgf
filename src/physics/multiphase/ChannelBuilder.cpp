#include "src/physics/multiphase/ChannelBuilder.H"
#include "src/physics/multiphase/MultiPhase.H"
#include "src/utilities/math_ops.H"
#include "src/utilities/IOManager.H"
#include "src/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_Gpu.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace kynema_sgf::channelbuilder {

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
        0.25_rt);
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
get_local_dimensions(
    const amrex::Real& x,
    const amrex::Real& y,
    const amrex::Real& z,
    const amrex::Real& start_x,
    const amrex::Real& start_y,
    const amrex::Real& start_z,
    const amrex::Real& end_x,
    const amrex::Real& end_y,
    const amrex::Real& end_z,
    const amrex::Real& dim0_s,
    const amrex::Real& dim1_s,
    const amrex::Real& dim2_s,
    const amrex::Real& dim0_e,
    const amrex::Real& dim1_e,
    const amrex::Real& dim2_e)

{
    // Interpolate dimensions based on position along segment
    const amrex::Real seg_length = std::sqrt(
        (end_x - start_x) * (end_x - start_x) +
        (end_y - start_y) * (end_y - start_y) +
        (end_z - start_z) * (end_z - start_z));
    const amrex::Real point_length = std::sqrt(
        (x - start_x) * (x - start_x) + (y - start_y) * (y - start_y) +
        (z - start_z) * (z - start_z));
    const amrex::Real t = point_length / seg_length;
    const amrex::Real dim0 = dim0_s + t * (dim0_e - dim0_s);
    const amrex::Real dim1 = dim1_s + t * (dim1_e - dim1_s);
    const amrex::Real dim2 = dim2_s + t * (dim2_e - dim2_s);

    return amrex::GpuArray<amrex::Real, 3>{{dim0, dim1, dim2}};
}

[[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::GpuArray<amrex::Real, 3>
transform_to_local_coordinates(
    const bool& do_translation,
    const bool& is_active,
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
    amrex::Real xp = x;
    amrex::Real yp = y;
    amrex::Real zp = z;
    if (do_translation) {
        xp -= start_x;
        yp -= start_y;
        zp -= start_z;
    }

    // Switch sign for active vs passive rotation
    const amrex::Real s = is_active ? 1.0_rt : -1.0_rt;

    // Rotate around z-axis based on xy component of direction
    const amrex::Real mag_xy = std::sqrt(a * a + b * b);
    const amrex::Real cos_theta_xy = a / mag_xy;
    const amrex::Real sin_theta_xy = b / mag_xy;

    const amrex::Real xpp = xp * cos_theta_xy - s * yp * sin_theta_xy;
    const amrex::Real ypp = s * xp * sin_theta_xy + yp * cos_theta_xy;

    // Rotate around y-axis based on z component
    const amrex::Real mag = std::sqrt(a * a + b * b + c * c);
    const amrex::Real cos_theta_xpz = mag_xy / mag;
    const amrex::Real sin_theta_xpz = c / mag;

    // Local coordinates: xloc along segment, yloc lateral, zloc vertical
    const amrex::Real xloc = xpp * cos_theta_xpz - s * zp * sin_theta_xpz;
    const amrex::Real zloc = s * xpp * sin_theta_xpz + zp * cos_theta_xpz;
    const amrex::Real yloc = ypp;

    return amrex::GpuArray<amrex::Real, 3>{{xloc, yloc, zloc}};
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

    amrex::Vector<ChannelSegmentType> h_type;
    amrex::Vector<ChannelVelocityProfile> h_velocity_profile;
    amrex::Vector<amrex::Real> h_dim0_s, h_dim1_s, h_dim2_s;
    amrex::Vector<amrex::Real> h_dim0_e, h_dim1_e, h_dim2_e;
    amrex::Vector<amrex::Array<amrex::Real, 3>> h_segment_start;
    amrex::Vector<amrex::Array<amrex::Real, 3>> h_segment_end;
    amrex::Vector<amrex::Real> h_flow_speed;

    for (const auto& lbl : labels) {
        const std::string key = identifier() + "." + lbl;
        amrex::ParmParse pp1(key);
        ChannelSegmentType type = ChannelSegmentType::Ellipse;
        amrex::Real dim0_s = 0.0_rt;
        amrex::Real dim1_s = 0.0_rt;
        amrex::Real dim2_s = 0.0_rt;
        amrex::Real dim0_e = 0.0_rt;
        amrex::Real dim1_e = 0.0_rt;
        amrex::Real dim2_e = 0.0_rt;
        amrex::Vector<amrex::Real> seg_start{0.0_rt, 0.0_rt, 0.0_rt};
        amrex::Vector<amrex::Real> seg_end{0.0_rt, 0.0_rt, 0.0_rt};

        std::string stype = "Ellipse";
        pp1.query("type", stype);
        if (stype == "Ellipse") {
            h_type.emplace_back(ChannelSegmentType::Ellipse);
            if (pp1.contains("diameter")) {
                pp1.get("diameter", dim0_s);
                dim1_s = dim0_s;
                dim0_e = dim0_s;
                dim1_e = dim0_s;
            } else if (pp1.contains("diameter_start")) {
                pp1.get("diameter_start", dim0_s);
                pp1.get("diameter_end", dim0_e);
                dim1_s = dim0_s;
                dim1_e = dim0_e;
            } else if (pp1.contains("horizontal_axis")) {
                pp1.get("horizontal_axis", dim0_s);
                pp1.get("vertical_axis", dim1_s);
                dim0_e = dim0_s;
                dim1_e = dim1_s;
            } else {
                pp1.get("horizontal_axis_start", dim0_s);
                pp1.get("vertical_axis_start", dim1_s);
                pp1.get("horizontal_axis_end", dim0_e);
                pp1.get("vertical_axis_end", dim1_e);
            }
        } else if (stype == "Trapezoid") {
            h_type.emplace_back(ChannelSegmentType::Trapezoid);
            if (pp1.contains("top_width")) {
                pp1.get("top_width", dim0_s);
                pp1.get("bottom_width", dim1_s);
                pp1.get("height", dim2_s);
            } else {
                pp1.get("top_width_start", dim0_s);
                pp1.get("bottom_width_start", dim1_s);
                pp1.get("height_start", dim2_s);
                pp1.get("top_width_end", dim0_e);
                pp1.get("bottom_width_end", dim1_e);
                pp1.get("height_end", dim2_e);
            }
        } else {
            amrex::Abort(
                "Invalid channel segment type specified: " + stype +
                ". Only 'Ellipse' and 'Trapezoid' are supported.");
        }
        pp1.getarr("segment_start_point", seg_start);
        pp1.getarr("segment_end_point", seg_end);

        h_dim0_s.emplace_back(dim0_s);
        h_dim1_s.emplace_back(dim1_s);
        h_dim2_s.emplace_back(dim2_s);
        h_dim0_e.emplace_back(dim0_e);
        h_dim1_e.emplace_back(dim1_e);
        h_dim2_e.emplace_back(dim2_e);
        h_segment_start.emplace_back(
            amrex::Array<amrex::Real, 3>{
                seg_start[0], seg_start[1], seg_start[2]});
        h_segment_end.emplace_back(
            amrex::Array<amrex::Real, 3>{seg_end[0], seg_end[1], seg_end[2]});

        stype = "Uniform";
        pp1.query("velocity_profile", stype);
        if (stype == "Uniform") {
            h_velocity_profile.emplace_back(ChannelVelocityProfile::Uniform);
        } else if (stype == "Linear") {
            h_velocity_profile.emplace_back(ChannelVelocityProfile::Linear);
        } else if (stype == "Parabolic") {
            h_velocity_profile.emplace_back(ChannelVelocityProfile::Parabolic);
        } else {
            amrex::Abort(
                "Invalid channel velocity profile specified: " + stype +
                ". Only 'Uniform', 'Linear', and 'Parabolic' are supported.");
        }
        amrex::Real flow_speed = 0.0_rt;
        pp1.get("flow_speed", flow_speed);
        h_flow_speed.emplace_back(flow_speed);
    }

    const int nseg = static_cast<int>(h_type.size());
    m_type.resize(nseg);
    m_dim0_s.resize(nseg);
    m_dim1_s.resize(nseg);
    m_dim2_s.resize(nseg);
    m_dim0_e.resize(nseg);
    m_dim1_e.resize(nseg);
    m_dim2_e.resize(nseg);
    m_segment_start.resize(nseg);
    m_segment_end.resize(nseg);
    m_velocity_profile.resize(nseg);
    m_flow_speed.resize(nseg);

    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, h_type.begin(), h_type.end(), m_type.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, h_dim0_s.begin(), h_dim0_s.end(),
        m_dim0_s.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, h_dim1_s.begin(), h_dim1_s.end(),
        m_dim1_s.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, h_dim2_s.begin(), h_dim2_s.end(),
        m_dim2_s.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, h_dim0_e.begin(), h_dim0_e.end(),
        m_dim0_e.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, h_dim1_e.begin(), h_dim1_e.end(),
        m_dim1_e.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, h_dim2_e.begin(), h_dim2_e.end(),
        m_dim2_e.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, h_segment_start.begin(),
        h_segment_start.end(), m_segment_start.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, h_segment_end.begin(), h_segment_end.end(),
        m_segment_end.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, h_velocity_profile.begin(),
        h_velocity_profile.end(), m_velocity_profile.begin());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, h_flow_speed.begin(), h_flow_speed.end(),
        m_flow_speed.begin());
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
    auto& velocity = m_repo.get_field("velocity");
    auto vel_arrs = velocity(level).arrays();
    auto& blank_mfab = m_terrain_blank(level);
    auto blank_arrs = blank_mfab.arrays();

    const int nseg = static_cast<int>(m_type.size());
    const ChannelSegmentType* type_ptr = m_type.data();
    const amrex::Real* dim0_s_ptr = m_dim0_s.data();
    const amrex::Real* dim1_s_ptr = m_dim1_s.data();
    const amrex::Real* dim2_s_ptr = m_dim2_s.data();
    const amrex::Real* dim0_e_ptr = m_dim0_e.data();
    const amrex::Real* dim1_e_ptr = m_dim1_e.data();
    const amrex::Real* dim2_e_ptr = m_dim2_e.data();
    const amrex::Array<amrex::Real, 3>* start_ptr = m_segment_start.data();
    const amrex::Array<amrex::Real, 3>* end_ptr = m_segment_end.data();
    const ChannelVelocityProfile* velocity_profile_ptr =
        m_velocity_profile.data();
    const amrex::Real* flow_speed_ptr = m_flow_speed.data();
    const bool multiphase = is_multiphase;
    const amrex::Real land_level = m_land_level;

    // Set all velocity to 0 for the sake of blanked cells
    velocity.setVal(0.0_rt);

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
            for (int seg = 0; seg < nseg; ++seg) {
                const auto& start = start_ptr[seg];
                const auto& end = end_ptr[seg];
                const auto& seg_type = type_ptr[seg];
                const auto& dim0_s = dim0_s_ptr[seg];
                const auto& dim1_s = dim1_s_ptr[seg];
                const auto& dim2_s = dim2_s_ptr[seg];
                const auto& dim0_e = dim0_e_ptr[seg];
                const auto& dim1_e = dim1_e_ptr[seg];
                const auto& dim2_e = dim2_e_ptr[seg];
                const auto& velocity_profile = velocity_profile_ptr[seg];
                const auto& flow_speed = flow_speed_ptr[seg];

                // Check if point is within bounding planes of segment start and
                // end
                if (is_point_within_planes(
                        x, y, z, start[0], start[1], start[2], end[0], end[1],
                        end[2])) {

                    // Get local dimensions at this point along segment
                    const auto local_dims = get_local_dimensions(
                        x, y, z, start[0], start[1], start[2], end[0], end[1],
                        end[2], dim0_s, dim1_s, dim2_s, dim0_e, dim1_e, dim2_e);
                    const amrex::Real dim0 = local_dims[0];
                    const amrex::Real dim1 = local_dims[1];
                    const amrex::Real dim2 = local_dims[2];

                    // Transform to local segment coordinates
                    const auto local_coords = transform_to_local_coordinates(
                        true, false, x, y, z, start[0], start[1], start[2],
                        end[0], end[1], end[2]);
                    const amrex::Real yloc = local_coords[1];
                    const amrex::Real zloc = local_coords[2];

                    // Transform flow speed to channel-aligned velocity vector
                    const auto local_flow = transform_to_local_coordinates(
                        false, true, flow_speed, 0.0_rt, 0.0_rt, start[0],
                        start[1], start[2], end[0], end[1], end[2]);
                    const amrex::Real uloc = local_flow[0];
                    const amrex::Real vloc = local_flow[1];
                    const amrex::Real wloc = local_flow[2];

                    if (seg_type == ChannelSegmentType::Ellipse) {
                        const bool inside_channel_segment =
                            ellipse(dim0, dim1, yloc, zloc);
                        outside_channel &= !inside_channel_segment;
                        // Set velocity according to profile if within channel
                        if (inside_channel_segment) {
                            amrex::Real speed_factor = 1.0_rt;
                            if (velocity_profile ==
                                ChannelVelocityProfile::Linear) {
                                const auto d =
                                    std::sqrt(yloc * yloc + zloc * zloc);
                                const auto d_ext = std::sqrt(
                                    utils::powi(
                                        (dim0 * yloc) / (2.0_rt * d), 2) +
                                    utils::powi(
                                        (dim1 * zloc) / (2.0_rt * d), 2));
                                speed_factor -= (d / d_ext);
                            } else if (
                                velocity_profile ==
                                ChannelVelocityProfile::Parabolic) {
                                const auto d =
                                    std::sqrt(yloc * yloc + zloc * zloc);
                                const auto d_ext = std::sqrt(
                                    utils::powi(
                                        (dim0 * yloc) / (2.0_rt * d), 2) +
                                    utils::powi(
                                        (dim1 * zloc) / (2.0_rt * d), 2));
                                speed_factor -= utils::powi(d / d_ext, 2);
                            } // Uniform case is default (else)
                            vel_arrs[nbx](i, j, k, 0) = uloc * speed_factor;
                            vel_arrs[nbx](i, j, k, 1) = vloc * speed_factor;
                            vel_arrs[nbx](i, j, k, 2) = wloc * speed_factor;
                        }
                    } else if (seg_type == ChannelSegmentType::Trapezoid) {
                        outside_channel &=
                            !trapezoid(dim0, dim1, dim2, yloc, zloc);
                    }
                }
            }

            // Do adjustments for multiphase case
            if (multiphase && z > land_level) {
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

} // namespace kynema_sgf::channelbuilder
