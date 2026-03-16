#ifndef CARTBOXREFINEMENT_H
#define CARTBOXREFINEMENT_H

#include <iosfwd>

#include "amr-wind/utilities/tagging/RefinementCriteria.H"
#include "AMReX_RealBox.H"
#include "AMReX_Vector.H"

namespace amr_wind {

/** Static refinement with Cartesian-aligned bounding boxes
 *  \ingroup amr_utils wind
 *
 *  Implements tagging functionality for nested refinement of meshes using
 *  bounding box specifications. This class only allows nested refinement of
 *  regions that are aligned with the coordinate directons.
 */
class CartBoxRefinement : public RefinementCriteria::Register<CartBoxRefinement>
{
public:
    static std::string identifier() { return "CartBoxRefinement"; }

    explicit CartBoxRefinement(CFDSim& sim);

    ~CartBoxRefinement() override = default;

    //! Read input file and initialize boxarray used to refine each level
    void initialize(const std::string& key) override;

    void
    operator()(int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

    //! Helper function to process inputs
    //!
    //! Created separate from initialize to allow unit testing
    void read_inputs(const amrex::AmrCore& /*mesh*/, std::istream& /*ifh*/);

    //! Vector of boxarrays that define refinement zones at each level
    [[nodiscard]] const amrex::Vector<amrex::BoxArray>& boxarray_vec() const
    {
        return m_boxarrays;
    }

protected:
    const amrex::AmrCore& m_mesh;

    //! Domain bounding boxes where refinement is performed at each level
    amrex::Vector<amrex::Vector<amrex::RealBox>> m_real_boxes;

    //! Boxarrays for each level in AMR hierarchy
    amrex::Vector<amrex::BoxArray> m_boxarrays;

    //! Number of levels of fixed nested refinement
    int m_nlevels{-1};
};

} // namespace amr_wind

#endif /* CARTBOXREFINEMENT_H */
