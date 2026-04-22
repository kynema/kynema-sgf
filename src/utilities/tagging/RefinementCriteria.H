#ifndef REFINEMENTCRITERIA_H
#define REFINEMENTCRITERIA_H

#include "AMReX_AmrCore.H"
#include "AMReX_TagBox.H"

#include "amr-wind/core/Factory.H"

/**
 *  \defgroup amr_utils Mesh refinement
 *  Mesh refinement
 *
 *  This group contains utilities used for static or adaptive refinement during
 *  a simulation.
 *
 *  \ingroup utilities
 */

namespace amr_wind {

class CFDSim;

/** Abstract interface for tagging cells for refinement
 *  \ingroup amr_utils
 *
 *  This class provides an API that can be used by concrete implementations to
 *  tag cells that must be refined based on a pre-defined criteria. The criteria
 *  can be a heuristic determined from the solution (e.g., vorticity magntiude
 *  or gradients), or some user defined criteria (e.g., static, nested
 * refinements).
 */
class RefinementCriteria : public Factory<RefinementCriteria, CFDSim&>
{
public:
    static std::string base_identifier() { return "RefinementCriteria"; }

    RefinementCriteria() = default;

    ~RefinementCriteria() override = default;

    //! Perform initialization tasks (e.g., read from file etc.). Called before
    //! mesh is initialized
    virtual void initialize(const std::string& key) = 0;

    /** Perform tagging operation based on refinement criteria
     *
     *  Interface is the same as amrex::AmrCore::ErrorEst
     */
    virtual void operator()(
        int level, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) = 0;
};

/** A collection of refinement criteria instances that are active during a
 * simulation \ingroup amr_utils
 */
class RefineCriteriaManager
{
public:
    explicit RefineCriteriaManager(CFDSim& sim);

    ~RefineCriteriaManager() = default;

    void initialize();

    void
    tag_cells(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow);

private:
    CFDSim& m_sim;

    amrex::Vector<std::unique_ptr<RefinementCriteria>> m_refiners;
};

} // namespace amr_wind

#endif /* REFINEMENTCRITERIA_H */
