#ifndef VORTICITYREFINEMENT_H
#define VORTICITYREFINEMENT_H

#include "amr-wind/utilities/tagging/RefinementCriteria.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {
class Field;

/** AMR refinement using vorticity
 *  \ingroup amr_utils
 *
 *  ```
 *  tagging.labels = v1
 *  tagging/v1.type = VorticityMagRefinement
 *  # Tagging based on value
 *  tagging/v1.values = 10.0_rt 10.0_rt

 *  ```
 */
class VorticityMagRefinement
    : public RefinementCriteria::Register<VorticityMagRefinement>
{
public:
    static std::string identifier() { return "VorticityMagRefinement"; }

    explicit VorticityMagRefinement(const CFDSim& sim);

    ~VorticityMagRefinement() override = default;

    //! Read input file and initialize boxarray used to refine each level
    void initialize(const std::string& key) override;

    void
    operator()(int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

private:
    const CFDSim& m_sim;

    Field* m_vel{nullptr};

    amrex::Vector<amrex::Real> m_vort_value;

    int m_max_lev_field{-1};
};

} // namespace amr_wind

#endif /* VORTICITYREFINEMENT_H */
