#ifndef CYLINDERREFINER_H
#define CYLINDERREFINER_H

#include "amr-wind/utilities/tagging/GeometryRefinement.H"
#include "amr-wind/core/vs/vector_space.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind::tagging {

/** Refine region inside a cylinder
 *  \ingroup amr_utils
 *
 *  The axis of the cylinder is defined by two vectors `start` and `end` that
 *  define the beginning and end of the cylinder along its axis. The outer
 *  extent is defined by `outer_radius`. An optional `inner_radius` can be
 *  specified to restrict the refinement to an annulus of the cylinder.
 */
class CylinderRefiner : public GeometryType::Register<CylinderRefiner>
{
public:
    static std::string identifier() { return "cylinder"; }

    CylinderRefiner(const CFDSim& /*unused*/, const std::string& key);

    ~CylinderRefiner() override = default;

    void operator()(
        const amrex::Box& /*bx*/,
        const amrex::Geometry& geom,
        const amrex::Array4<amrex::TagBox::TagType>& tags) const override;
    [[nodiscard]] const amrex::RealBox& bound_box() const override
    {
        return m_bound_box;
    };

private:
    vs::Vector m_start;
    vs::Vector m_end;

    amrex::Real m_outer_radius;
    amrex::Real m_inner_radius{0.0_rt};

    amrex::RealBox m_bound_box;
};

} // namespace amr_wind::tagging

#endif /* CYLINDERREFINER_H */
