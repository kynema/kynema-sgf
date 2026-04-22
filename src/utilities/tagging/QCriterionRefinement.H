#ifndef QCRITERIONREFINEMENT_H
#define QCRITERIONREFINEMENT_H

#include "amr-wind/utilities/tagging/RefinementCriteria.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {
class Field;

/** AMR refinement using Q criterion
 *  \ingroup amr_utils
 *
 *  ```
 *  tagging.labels = qc1
 *  tagging/qc1.type = QCriterionRefinement
 *  # Tagging based on value
 *  tagging/qc1.values = 10.0_rt 10.0_rt

 *  ```
 */
class QCriterionRefinement
    : public RefinementCriteria::Register<QCriterionRefinement>
{
public:
    static std::string identifier() { return "QCriterionRefinement"; }

    explicit QCriterionRefinement(const CFDSim& sim);

    ~QCriterionRefinement() override = default;

    //! Read input file and initialize boxarray used to refine each level
    void initialize(const std::string& key) override;

    void
    operator()(int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

private:
    const CFDSim& m_sim;

    Field* m_vel{nullptr};

    amrex::Vector<amrex::Real> m_qc_value;

    int m_max_lev_field{-1};

    bool m_nondim{true};
};

} // namespace amr_wind

#endif /* QCRITERIONREFINEMENT_H */
