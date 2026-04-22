#ifndef OVERSETOPS_H
#define OVERSETOPS_H

#include "amr-wind/CFDSim.H"
#include "amr-wind/physics/multiphase/MultiPhase.H"
#include "AMReX_REAL.H"

using namespace amrex::literals;

namespace amr_wind {

class OversetOps
{
public:
    OversetOps() = default;

    void initialize(CFDSim& sim);
    void pre_advance_work();
    void post_advance_work();

    void update_gradp();

private:
    // Functions called within public functions
    void parameter_output() const;
    void form_perturb_pressure();
    void replace_masked_gradp();

    // Check for multiphase sim
    bool m_vof_exists{false};

    // Coupling options
    bool m_replace_gradp_postsolve{true};

    // Verbosity
    int m_verbose{0};

    // Pointer for pressure gradient copy field
    amr_wind::Field* m_gp_copy{nullptr};
    // Pointer for MultiPhase physics
    amr_wind::MultiPhase* m_mphase{nullptr};

    CFDSim* m_sim_ptr;
};

} // namespace amr_wind

#endif /* OVERSETOPS_H */
