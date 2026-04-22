#ifndef BOXREFINER_H
#define BOXREFINER_H

#include "amr-wind/utilities/tagging/GeometryRefinement.H"
#include "amr-wind/core/vs/vector_space.H"

namespace amr_wind::tagging {

/** AMR inside an arbitrarily oriented hexahedral box
 *  \ingroup amr_utils
 *
 *  The hexahedral box is defined by three vectors `xaxis`, `yaxis`, and `zaxis`
 *  that specify the direction and lengths of the box with respect to `origin`.
 *  If the axes are aligned with Cartesian coordinates, it is better to use the
 *  CartBoxRefinement class for faster tagging.
 *
 *  \sa CartBoxRefinement
 */
class BoxRefiner : public GeometryType::Register<BoxRefiner>
{
public:
    static std::string identifier() { return "box"; }

    BoxRefiner(const CFDSim& /*unused*/, const std::string& key);

    ~BoxRefiner() override = default;

    void operator()(
        const amrex::Box& /*bx*/,
        const amrex::Geometry& geom,
        const amrex::Array4<amrex::TagBox::TagType>& tags) const override;

    [[nodiscard]] const amrex::RealBox& bound_box() const override
    {
        return m_bound_box;
    };

protected:
    amrex::Gpu::DeviceVector<vs::Vector> m_hex_corners;
    amrex::Gpu::DeviceVector<vs::Vector> m_face_normals;
    amrex::Gpu::DeviceVector<int> m_face_origin;
    amrex::RealBox m_bound_box;
};

} // namespace amr_wind::tagging

#endif /* BOXREFINER_H */
