#ifndef OVERSETREFINEMENT_H
#define OVERSETREFINEMENT_H

#include "src/utilities/tagging/RefinementCriteria.H"

namespace kynema_sgf {

class OversetRefinement : public RefinementCriteria::Register<OversetRefinement>
{
public:
    static std::string identifier() { return "OversetRefinement"; }

    explicit OversetRefinement(const CFDSim& sim);

    ~OversetRefinement() override = default;

    //! Read input file and initialize boxarray used to refine each level
    void initialize(const std::string& key) override;

    void
    operator()(int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
        override;

private:
    const CFDSim& m_sim;

    int m_max_lev;
    bool m_tag_fringe{false};
    bool m_tag_hole{false};
};

} // namespace kynema_sgf

#endif /* OVERSETREFINEMENT_H */
