#ifndef TIOGAINTERFACE_H
#define TIOGAINTERFACE_H

#include <vector>
#include "amr-wind/overset/OversetManager.H"
#include "amr-wind/overset/overset_types.H"

namespace amr_wind {

class IntField;
class ScratchField;
class IntScratchField;

/** AMR mesh data in TIOGA format
 */
struct AMROversetInfo
{
    template <typename T>
    using AType = tioga_iface::AmrDualArray<T>;

    /**
     *  \param nglobal Total number of patches
     *  \param nlocal  Number of patches on the current MPI rank
     */
    AMROversetInfo(int nglobal, int nlocal);

    // Arrays of size ngrids_global
    AType<int> level;
    AType<int> mpi_rank;
    AType<int> local_id;
    AType<int> ilow;
    AType<int> ihigh;
    AType<int> dims;
    AType<amrex::Real> xlo;
    AType<amrex::Real> dx;

    // Arrays of size ngrids_local
    AType<int> global_idmap;
    AType<int*> iblank_node;
    AType<int*> iblank_cell;
    AType<amrex::Real*> qcell;
    AType<amrex::Real*> qnode;

    const int ngrids_global;
    const int ngrids_local;
};

class TiogaInterface : public OversetManager::Register<TiogaInterface>
{
public:
    static std::string identifier() { return "TIOGA"; }

    explicit TiogaInterface(CFDSim& sim);

    /** Perform one-time initialization actions.
     *
     *  This method is called after the initial mesh and solution fields have
     *  been created.
     */
    void post_init_actions() override;

    /** Perform necessary updates after a regrid.
     */
    void post_regrid_actions() override;

    /** Perform actions before initializing connectivity through overset
     * library.
     *
     *  This method is called before performing connectivity. The implementation
     *  should update all necessary mesh data and metadata that will be used for
     *  overset connectivity.
     */
    void pre_overset_conn_work() override;

    /** Updates necessary after overset connectivity has been updated.
     *
     *  Update dependent fields or data structures as a result of changes in
     *  overset connectivity.
     */
    void post_overset_conn_work() override;

    /** Perform actions necessary to register fields for overset exchange.
     *
     *  This method will be called at least once every timestep and, perhaps,
     *  multiple times within a timestep.
     */
    void register_solution(
        const std::vector<std::string>& cell_vars,
        const std::vector<std::string>& node_vars) override;

    /** Additional work after solution fields have been exchanged.
     */
    void update_solution() override;

    AMROversetInfo& amr_overset_info() { return *m_amr_data; }

    ScratchField& qvars_cell()
    {
        AMREX_ALWAYS_ASSERT(m_qcell_host);
        return *m_qcell_host;
    }
    ScratchField& qvars_node()
    {
        AMREX_ALWAYS_ASSERT(m_qnode_host);
        return *m_qnode_host;
    }

private:
    void amr_to_tioga_mesh();

    void amr_to_tioga_iblank();

    CFDSim& m_sim;

    //! IBLANK on cell centered fields
    IntField& m_iblank_cell;

    //! IBLANK for nodal fields
    IntField& m_iblank_node;

    //! IBLANK on cell centered fields on host
    std::unique_ptr<IntScratchField> m_iblank_cell_host;

    //! IBLANK for nodal fields on host
    std::unique_ptr<IntScratchField> m_iblank_node_host;

    //! AMReX mask for linear system solves
    IntField& m_mask_cell;

    IntField& m_mask_node;

    std::unique_ptr<AMROversetInfo> m_amr_data;

    std::unique_ptr<ScratchField> m_qcell;
    std::unique_ptr<ScratchField> m_qnode;

    std::unique_ptr<ScratchField> m_qcell_host;
    std::unique_ptr<ScratchField> m_qnode_host;

    std::vector<std::string> m_cell_vars;
    std::vector<std::string> m_node_vars;
};

} // namespace amr_wind

#endif /* TIOGAINTERFACE_H */
