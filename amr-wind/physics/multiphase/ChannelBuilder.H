#ifndef CHANNELBUILDER_H
#define CHANNELBUILDER_H

#include "amr-wind/core/Physics.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/CFDSim.H"

enum class ChannelSegmentType : std::uint8_t { Ellipse, Trapezoid };

namespace amr_wind::channelbuilder {

[[nodiscard]] bool trapezoid(
    const amrex::Real& top,
    const amrex::Real& bottom,
    const amrex::Real& height,
    const amrex::Real& hcoord,
    const amrex::Real& vcoord);

[[nodiscard]] bool ellipse(
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

/** Transform global coordinates to local segment coordinates.
 *  Returns [xloc, yloc, zloc] where xloc is along segment,
 *  yloc is lateral/horizontal, zloc is vertical to segment.
 */
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
    const amrex::Real& end_z);

/** Terraindrag Flow physics
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

    amrex::Real m_water_level{0.0_rt};
    amrex::Real m_land_level{0.0_rt};

    //! Quantities for constructing channels
    amrex::Vector<ChannelSegmentType> m_type;
    amrex::Vector<amrex::Real> m_dim0;
    amrex::Vector<amrex::Real> m_dim1;
    amrex::Vector<amrex::Real> m_dim2;
    amrex::Vector<amrex::Array<amrex::Real, 3>> m_segment_start;
    amrex::Vector<amrex::Array<amrex::Real, 3>> m_segment_end;

    bool is_multiphase;
};
} // namespace amr_wind::channelbuilder

#endif
