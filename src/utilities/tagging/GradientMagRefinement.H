#ifndef GRADIENTMAGREFINEMENT_H
#define GRADIENTMAGREFINEMENT_H

#include "amr-wind/utilities/tagging/RefinementCriteria.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {
class Field;

/** AMR refinement using a given field (e.g., density)
 *  \ingroup amr_utils
 *
 *  ```
 *  tagging.labels = t1
 *  tagging/t1.type = GradientMagRefinement
 *  tagging/t1.field_name = density
 *  # Tagging based on field's gradient magnitude value
 *  # tagging/t1.values = 10.0_rt 10.0_rt
 *  ```
 */
class GradientMagRefinement
    : public RefinementCriteria::Register<GradientMagRefinement>
{
public:
    static std::string identifier() { return "GradientMagRefinement"; }

    explicit GradientMagRefinement(const CFDSim& sim);

    ~GradientMagRefinement() override = default;

    //! Read input file and initialize boxarray used to refine each level
    void initialize(const std::string& key) override;

    void
    operator()(int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

private:
    const CFDSim& m_sim;

    Field* m_field{nullptr};

    amrex::Vector<amrex::Real> m_gradmag_value;

    int m_max_lev_field{-1};
};

} // namespace amr_wind

#endif /* GRADIENTMAGREFINEMENT_H */
