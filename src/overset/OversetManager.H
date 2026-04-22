#ifndef OVERSETMANAGER_H
#define OVERSETMANAGER_H

#include "amr-wind/core/Factory.H"

namespace amr_wind {

class CFDSim;

/**
 *  \defgroup overset Overset connectivity interface
 *
 *  Overset connectivity interface.
 *
 *  This module contains all the code necessary for performing overset
 *  simulations by integrating AMR-Wind with other solvers using overset meshes.
 */

/** Overset interface manager
 *  \ingroup overset
 */
class OversetManager : public Factory<OversetManager, CFDSim&>
{
public:
    static std::string base_identifier() { return "OversetManager"; }

    /** Perform one-time initialization actions.
     *
     *  This method is called after the initial mesh and solution fields have
     *  been created.
     */
    virtual void post_init_actions() = 0;

    /** Perform necessary updates after a regrid.
     */
    virtual void post_regrid_actions() = 0;

    /** Perform actions before initializing connectivity through overset
     * library.
     *
     *  This method is called before performing connectivity. The implementation
     *  should update all necessary mesh data and metadata that will be used for
     *  overset connectivity.
     */
    virtual void pre_overset_conn_work() = 0;

    /** Updates necessary after overset connectivity has been updated.
     *
     *  Update dependent fields or data structures as a result of changes in
     *  overset connectivity.
     */
    virtual void post_overset_conn_work() = 0;

    /** Perform actions necessary to register fields for overset exchange.
     *
     *  This method will be called at least once every timestep and, perhaps,
     *  multiple times within a timestep.
     */
    virtual void register_solution(
        const std::vector<std::string>& cell_vars,
        const std::vector<std::string>& node_vars) = 0;

    /** Additional work after solution fields have been exchanged.
     */
    virtual void update_solution() = 0;
};

} // namespace amr_wind

#endif /* OVERSETMANAGER_H */
