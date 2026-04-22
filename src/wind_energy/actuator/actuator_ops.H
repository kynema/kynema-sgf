#ifndef ACT_OPS_H
#define ACT_OPS_H

#include "amr-wind/wind_energy/actuator/actuator_types.H"
#include "amr-wind/wind_energy/actuator/ActParser.H"
#include "AMReX_Vector.H"

namespace amr_wind::actuator::ops {

/** Read user inputs from an input file for a given actuator object.
 *
 *  \ingroup actuator
 *
 *  This operator is called within Physics::pre_init_actions() and, therefore,
 *  does not have access to mesh data. Actions that require mesh information for
 *  initialization should be deferred until ops::init_data_structures call.
 */
template <typename ActTrait, typename SrcTrait, typename = void>
struct ReadInputsOp;

/** Perform one-time initialization of the actuator object.
 *
 *  \ingroup actuator
 *
 *  This operator is responsible for allocating memory ensuring that the
 *  actuator instance has been initialized and ready to be used by the Actuator
 *  physics instance. This function is called once during
 * Actuator::post_init_actions
 *
 */
template <typename ActTrait, typename SrcTrait, typename = void>
struct InitDataOp;

/** Update the positions of the actuator nodes.
 *
 *  \ingroup actuator
 *
 *  After this call, the interface expects updated position vectors in the
 *  vel_pos array of the actuator grid instance. These positions will be used
 *  the sample the velocity field during the simulation.
 */
template <typename ActTrait, typename SrcTrait, typename = void>
struct UpdatePosOp;

/** Process updated velocities from CFD flow field.
 *
 *  \ingroup actuator
 *
 *  The actuator physics interface will populate the `vel` array of the
 *  actuator grid instance before calling this function. The updated velocities
 *  can be further processed by specializing this function.
 */
template <typename ActTrait, typename SrcTrait, typename = void>
struct UpdateVelOp;

/** Compute aerodynamic forces at the actuator grid points during a simulation.
 *
 *  \ingroup actuator
 */
template <typename ActTrait, typename SrcTrait, typename = void>
struct ComputeForceOp;

/** Compute source term for the momentum equation.
 *
 *  \ingroup actuator
 */
template <typename ActTrait, typename SrcTrait, typename = void>
class ActSrcOp;

template <typename ActTrait, typename SrcTrait, typename = void>
struct ProcessOutputsOp;

/** Determine the list of processes (MPI ranks) that is expected to be
 *  influenced by this actuator source.
 *
 *  The list of processes is used to determine where velocity sampling data is
 *  too be gathered, as well as ensuring that all actuator data is properly
 *  mapped.
 */
template <typename T>
void determine_influenced_procs(typename T::DataType& /*data*/);

/** During initialization determine the influenced processes and elect a root
 *  process to manage I/O and external solver interfaces for this turbine.
 *
 *  \tparam T An actuator traits type
 *  \param  data Data object for the specific actuator instance
 *
 *  \param act_proc_count A map of the number of turbines managed by each proc
 */
template <typename T>
void determine_root_proc(
    typename T::DataType& /*data*/, amrex::Vector<int>& /*act_proc_count*/);

} // namespace amr_wind::actuator::ops

#include "amr-wind/wind_energy/actuator/actuator_opsI.H"

#endif /* ACT_OPS_H */
