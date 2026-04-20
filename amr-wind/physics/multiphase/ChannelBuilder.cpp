#include "amr-wind/physics/multiphase/ChannelBuilder.H"
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

namespace amr_wind::channelbuilder {

namespace {} // namespace

ChannelBuilder::ChannelBuilder(CFDSim& sim)
    : m_sim(sim)
    , m_repo(sim.repo())
    , m_mesh(sim.mesh())
    , m_terrain_blank(sim.repo().declare_int_field("terrain_blank", 1, 1, 1))
{

    m_sim.io_manager().register_output_int_var("terrain_drag");

    m_terrain_blank.setVal(0);
    amrex::Vector<std::string> labels;
    amrex::ParmParse pp(identifier());
    is_multiphase = pp.contains("water_level");
    if (is_multiphase) {
        pp.get("water_level", m_water_level);
        pp.get("land_level", m_land_level);
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
            } else {
                pp1.get("lateral_axis", dim0);
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

        m_type.emplace_back(stype);
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
