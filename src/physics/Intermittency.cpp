#include "src/physics/Intermittency.H"
#include "src/CFDSim.H"
#include "AMReX_ParmParse.H"
#include "src/utilities/trig_ops.H"

namespace amr_wind {

Intermittency::Intermittency(const CFDSim& sim)
    : m_intermittency(sim.repo().declare_field(
          "intermittency", 1, sim.pde_manager().num_ghost_state()))
{
    amrex::ParmParse pp("incflo");
    pp.query("gamma_intermittency", m_gamma);
    auto& gamma_int_fld = sim.repo().get_field("intermittency");
    gamma_int_fld.set_default_fillpatch_bc(sim.time());
}

/** Initialize the intermittency field at the beginning of the
 *  simulation.
 */
void Intermittency::initialize_fields(
    int level, const amrex::Geometry& /*geom*/)
{
    auto& gamma = m_intermittency(level);

    gamma.setVal(m_gamma);
}

} // namespace amr_wind
