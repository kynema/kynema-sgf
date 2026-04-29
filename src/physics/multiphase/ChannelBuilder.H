#ifndef CHANNELBUILDER_H
#define CHANNELBUILDER_H

#include "src/core/Physics.H"
#include "src/core/Field.H"
#include "src/CFDSim.H"
#include "AMReX_Gpu.H"

namespace kynema_sgf::channelbuilder {

enum class ChannelSegmentType : std::uint8_t { Ellipse, Trapezoid };

enum class ChannelVelocityProfile : std::uint8_t { Uniform, Linear, Parabolic };
[[nodiscard]] AMREX_GPU_HOST_DEVICE bool trapezoid(
    const amrex::Real& top,
    const amrex::Real& bottom,
    const amrex::Real& height,
    const amrex::Real& hcoord,
    const amrex::Real& vcoord);

[[nodiscard]] AMREX_GPU_HOST_DEVICE bool ellipse(
    const amrex::Real& ax_horz,
    const amrex::Real& ax_vert,
    const amrex::Real& hcoord,
    const amrex::Real& vcoord);

/** Check if a point is within the bounding planes of a segment
 *  defined by start and end points.
 */
[[nodiscard]] AMREX_GPU_HOST_DEVICE bool is_point_within_planes(
    const amrex::Real& x,
    const amrex::Real& y,
    const amrex::Real& z,
    const amrex::Real& start_x,
    const amrex::Real& start_y,
    const amrex::Real& start_z,
    const amrex::Real& end_x,
    const amrex::Real& end_y,
    const amrex::Real& end_z);

/** Calculate local shape dimensions according to linear interpolation between
 * start and end dimensions using length along segment.
 */
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
    const amrex::Real& dim2_e);

/** Transform global coordinates to local segment coordinates.
 *  Returns [xloc, yloc, zloc] where xloc is along segment,
 *  yloc is lateral/horizontal, zloc is vertical to segment.
 */
[[nodiscard]] AMREX_GPU_HOST_DEVICE amrex::GpuArray<amrex::Real, 3>
transform_to_local_coordinates(
    bool do_translation,
    bool is_active,
    const amrex::Real& x,
    const amrex::Real& y,
    const amrex::Real& z,
    const amrex::Real& start_x,
    const amrex::Real& start_y,
    const amrex::Real& start_z,
    const amrex::Real& end_x,
    const amrex::Real& end_y,
    const amrex::Real& end_z);

/** Channel segment geometry and velocity initialization physics.
 *  \ingroup physics
 */

class ChannelBuilder : public Physics::Register<ChannelBuilder>
{
public:
    static std::string identifier() { return "ChannelBuilder"; }

    explicit ChannelBuilder(CFDSim& sim);

    ~ChannelBuilder() override = default;

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

    //! Blanking field for terrain or buildings
    IntField& m_terrain_blank;

    bool m_is_multiphase{false};
    amrex::Real m_rho_init{1.0_rt};
    amrex::Real m_water_level{0.0_rt};
    amrex::Real m_land_level{0.0_rt};

    //! Quantities for constructing channels
    amrex::Gpu::DeviceVector<ChannelSegmentType> m_type;
    amrex::Gpu::DeviceVector<amrex::Real> m_dim0_s;
    amrex::Gpu::DeviceVector<amrex::Real> m_dim1_s;
    amrex::Gpu::DeviceVector<amrex::Real> m_dim2_s;
    amrex::Gpu::DeviceVector<amrex::Real> m_dim0_e;
    amrex::Gpu::DeviceVector<amrex::Real> m_dim1_e;
    amrex::Gpu::DeviceVector<amrex::Real> m_dim2_e;
    amrex::Gpu::DeviceVector<amrex::Array<amrex::Real, 3>> m_segment_start;
    amrex::Gpu::DeviceVector<amrex::Array<amrex::Real, 3>> m_segment_end;

    //! Quantities for setting velocity in channels
    amrex::Gpu::DeviceVector<ChannelVelocityProfile> m_velocity_profile;
    amrex::Gpu::DeviceVector<amrex::Real> m_flow_speed;
};
} // namespace kynema_sgf::channelbuilder

#endif
