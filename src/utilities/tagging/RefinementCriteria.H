#ifndef REFINEMENTCRITERIA_H
#define REFINEMENTCRITERIA_H

#include "AMReX_AmrCore.H"
#include "AMReX_TagBox.H"

#include "src/core/Factory.H"

/**
 *  \defgroup amr_utils Mesh refinement
 *  Mesh refinement
 *
 *  This group contains utilities used for static or adaptive refinement during
 *  a simulation.
 *
 *  \ingroup utilities
 */

namespace kynema_sgf {

namespace tagging {

/** Enumeration for tagging operation types */
enum class TaggingOperator : std::uint8_t { AND, OR, AND_NOT, OR_NOT };

/** Convert string to TaggingOperator enum
 *  \param op The operator string (case-sensitive), needs lowercase
 *  \return The corresponding TaggingOperator enum value
 *  \throw Calls amrex::Abort if the string is invalid
 */
AMREX_FORCE_INLINE TaggingOperator string_to_operator(const std::string& op)
{
    if (op == "and") {
        return TaggingOperator::AND;
    }
    if (op == "or") {
        return TaggingOperator::OR;
    }
    if (op == "and_not") {
        return TaggingOperator::AND_NOT;
    }
    if (op == "or_not") {
        return TaggingOperator::OR_NOT;
    }
    amrex::Abort(
        "Invalid operator for tagging criteria: " + op +
        ". Must be one of: and, or, and_not, or_not.");
    // Unreachable, but avoids compiler warning
    return TaggingOperator::AND;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::TagBox::TagVal tag_val(
    const bool previous_val, const bool current_val, const TaggingOperator& op)
{
    switch (op) {
    case TaggingOperator::AND:
        return (previous_val && current_val) ? amrex::TagBox::SET
                                             : amrex::TagBox::CLEAR;
    case TaggingOperator::OR:
        return (previous_val || current_val) ? amrex::TagBox::SET
                                             : amrex::TagBox::CLEAR;
    case TaggingOperator::AND_NOT:
        return (previous_val && !current_val) ? amrex::TagBox::SET
                                              : amrex::TagBox::CLEAR;
    case TaggingOperator::OR_NOT:
        return (previous_val || !current_val) ? amrex::TagBox::SET
                                              : amrex::TagBox::CLEAR;
    }
    // Unreachable, but avoids compiler warning
    return amrex::TagBox::CLEAR;
}

} // namespace tagging

class CFDSim;

/** Abstract interface for tagging cells for refinement
 *  \ingroup amr_utils
 *
 *  This class provides an API that can be used by concrete implementations to
 *  tag cells that must be refined based on a pre-defined criteria. The criteria
 *  can be a heuristic determined from the solution (e.g., vorticity magnitude
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

    //! Inactive criteria are skipped rather than combined as false predicates.
    [[nodiscard]] virtual bool active(int /*level*/) const { return true; }

    /** Accumulate matching cells into a predicate mask using TagBox::SET.
     *
     *  The manager supplies a cleared mask and combines it with existing tags.
     *  Multiple matches within one criterion form a union.
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

    struct RefinementStep
    {
        std::unique_ptr<RefinementCriteria> criterion;
        tagging::TaggingOperator op;
    };

    amrex::Vector<RefinementStep> m_refiners;
};

} // namespace kynema_sgf

#endif /* REFINEMENTCRITERIA_H */
