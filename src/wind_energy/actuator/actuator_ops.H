#ifndef ACT_OPS_H
#define ACT_OPS_H

#include "src/wind_energy/actuator/actuator_types.H"
#include "src/wind_energy/actuator/ActParser.H"
#include "AMReX_Vector.H"

namespace kynema_sgf::actuator::ops {

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
 *  Actuator::post_init_actions().
 *
 */
template <typename ActTrait, typename SrcTrait, typename = void>
struct InitDataOp;

/** Update the positions of the actuator nodes.
 *
 *  \ingroup actuator
 *
 *  After this call, the actuator grid's vel_pos array contains the updated
 *  position vectors. The Actuator physics-level method coordinates this
 *  actuator-type operation and updates the sampling container.
 */
template <typename ActTrait, typename SrcTrait, typename = void>
struct UpdatePosOp;

/** Update actuator fields from sampled flow data.
 *
 *  \ingroup actuator
 *
 *  ActuatorModel::update_fields() copies sampled velocity and density into the
 *  actuator grid, then dispatches to this actuator-type operation to update
 *  instance data.
 */
template <typename ActTrait, typename SrcTrait, typename = void>
struct UpdateFieldsOp;

/** Update the model velocity if it is prescribed (as opposed to being solved)
 *
 *  \ingroup actuator
 */
template <typename ActTrait, typename SrcTrait, typename = void>
struct PrescribeModelVelocityOp;

/** Compute actuator loads and scatter them to the relevant processes.
 *
 *  \ingroup actuator
 */
template <typename ActTrait, typename SrcTrait, typename = void>
struct ComputeLoadsOp;

/** Advance the actuator model by one CFD time step then scatter new position.
 *
 *  \ingroup actuator
 */
template <typename ActTrait, typename SrcTrait, typename = void>
struct AdvanceModelOp;

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
 *  The list of processes determines where sampled fields are gathered and
 *  ensures that all actuator data is properly mapped.
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

} // namespace kynema_sgf::actuator::ops

#include "src/wind_energy/actuator/actuator_opsI.H"

#endif /* ACT_OPS_H */
