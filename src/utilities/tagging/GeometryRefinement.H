#ifndef GEOMETRYREFINEMENT_H
#define GEOMETRYREFINEMENT_H

#include "amr-wind/utilities/index_operations.H"
#include "amr-wind/utilities/tagging/RefinementCriteria.H"

namespace amr_wind {
namespace tagging {

/** Abstract definition of shapes for AMR
 *  \ingroup amr_utils
 */
class GeometryType
    : public Factory<GeometryType, const CFDSim&, const std::string&>
{
public:
    static std::string base_identifier() { return "GeometryType"; }

    GeometryType() = default;

    ~GeometryType() override = default;

    virtual void operator()(
        const amrex::Box&,
        const amrex::Geometry& geom,
        const amrex::Array4<amrex::TagBox::TagType>& tags) const = 0;

    [[nodiscard]] virtual const amrex::RealBox& bound_box() const = 0;
};

} // namespace tagging

/** AMR based on geometric shapes
 *  \ingroup amr_utils
 */
class GeometryRefinement
    : public RefinementCriteria::Register<GeometryRefinement>
{
public:
    static std::string identifier() { return "GeometryRefinement"; }

    explicit GeometryRefinement(const CFDSim& /*sim*/);

    ~GeometryRefinement() override = default;

    //! Read input file and initialize boxarray used to refine each level
    void initialize(const std::string& key) override;

    void
    operator()(int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

private:
    const CFDSim& m_sim;

    amrex::Vector<std::unique_ptr<tagging::GeometryType>> m_geom_refiners;

    /** Option to refine a specific level
     *  If set greater than -1, then only that level is acted on
     */
    int m_set_level{-1};

    // If the user hasn't specified a particular level, act on a range of levels
    int m_min_level{0};
    int m_max_level{32};
};

} // namespace amr_wind

#endif /* GEOMETRYREFINEMENT_H */
