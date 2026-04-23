#ifndef CFDSIM_H
#define CFDSIM_H

#include "AMReX_AmrCore.H"
#include "src/core/SimTime.H"
#include "src/core/FieldRepo.H"
#include "src/equation_systems/PDEBase.H"
#include "src/core/Physics.H"
#include "src/boundary_conditions/field_boundary_fill/FieldBoundary.H"
#include "src/core/MeshMap.H"
#include "src/helics.H"

#include "src/boundary_conditions/field_boundary_fill/ReadERFFunction.H"
class MultiBlockContainer;

/** Kynema-SGF
 *
 *  All C++ code in Kynema-SGF is organized within the kynema_sgf namespace.
 */
namespace kynema_sgf {
class IOManager;
class PostProcessManager;
class OversetManager;
class ExtSolverMgr;
class HelicsStorage;

namespace transport {
class TransportModel;
}

namespace turbulence {
class TurbulenceModel;
}

/**
 *  \defgroup core Kynema-SGF core
 *  Core data structures in Kynema-SGF
 *
 *  This group describes the core data structures implemented in Kynema-SGF that
 *  provide a higher-level abstraction over the classes available in AMReX.
 */

/** Data structures for a CFD simulation
 *  \ingroup core
 *
 *  CFDSim represents the kynema-sgf execution environment and manages all the
 *  necessary components used to perform a simulation. Each simulation contains
 *  a single CFDSim instance that holds references to the mesh, FieldRepo,
 *  SimTime, PhysicsMgr, pde::PDEMgr, IOManager, and post-processing manager
 *  instances. This class is just a data holder and does not perform any
 *  computational logic itself. The computational heavy-lifting is performed by
 *  the class instances within a time integration loop.
 */
class CFDSim
{
public:
    explicit CFDSim(amrex::AmrCore& mesh);

    ~CFDSim();

    //! Return the AMR mesh hierarchy
    amrex::AmrCore& mesh() { return m_mesh; }
    const amrex::AmrCore& mesh() const { return m_mesh; }

    //! Return simulation time control
    SimTime& time() { return m_time; }
    const SimTime& time() const { return m_time; }

    void set_mbc(MultiBlockContainer* mbc) { m_mbc = mbc; }
    MultiBlockContainer** mbc() { return &m_mbc; }

    void set_read_erf(ReadERFFunction fcn) { m_read_erf = std::move(fcn); }
    ReadERFFunction* get_read_erf() { return &m_read_erf; }

    //! Return the field repository
    FieldRepo& repo() { return m_repo; }
    FieldRepo& repo() const { return m_repo; }

    pde::PDEMgr& pde_manager() { return m_pde_mgr; }
    const pde::PDEMgr& pde_manager() const { return m_pde_mgr; }

    PhysicsMgr& physics_manager() { return m_physics_mgr; }
    const PhysicsMgr& physics_manager() const { return m_physics_mgr; }

    PhysicsMgr::TypeVector& physics() { return m_physics_mgr.objects(); }
    const PhysicsMgr::TypeVector& physics() const
    {
        return m_physics_mgr.objects();
    }

    FieldBoundaryMgr& field_boundary_manager() { return m_field_boundary_mgr; }
    const FieldBoundaryMgr& field_boundary_manager() const
    {
        return m_field_boundary_mgr;
    }

    FieldBoundaryMgr::TypeVector& field_boundaries()
    {
        return m_field_boundary_mgr.objects();
    }
    const FieldBoundaryMgr::TypeVector& field_boundaries() const
    {
        return m_field_boundary_mgr.objects();
    }

    transport::TransportModel& transport_model() { return *m_transport; }
    const transport::TransportModel& transport_model() const
    {
        return *m_transport;
    }

    turbulence::TurbulenceModel& turbulence_model() { return *m_turbulence; }
    const turbulence::TurbulenceModel& turbulence_model() const
    {
        return *m_turbulence;
    }

    IOManager& io_manager() { return *m_io_mgr; }
    const IOManager& io_manager() const { return *m_io_mgr; }

    PostProcessManager& post_manager() { return *m_post_mgr; }
    const PostProcessManager& post_manager() const { return *m_post_mgr; }

    OversetManager* overset_manager() { return m_overset_mgr.get(); }
    const OversetManager* overset_manager() const
    {
        return m_overset_mgr.get();
    }

    MeshMap* mesh_mapping() { return m_mesh_map.get(); }
    const MeshMap* mesh_mapping() const { return m_mesh_map.get(); }

    ExtSolverMgr& ext_solver_manager() { return *m_ext_solver_mgr; }
    const ExtSolverMgr& ext_solver_manager() const { return *m_ext_solver_mgr; }

    HelicsStorage& helics() { return *m_helics; }
    const HelicsStorage& helics() const { return *m_helics; }

    bool has_overset() const;

    //! Instantiate the transport model based on user inputs
    void create_transport_model();

    //! Instantiate the turbulence model based on user inputs
    void create_turbulence_model();

    //! Instantiate all physics instances based on user inputs
    void init_physics();

    //! Instantiate all field boundaries based on user inputs
    void init_field_boundaries();

    //! Activate overset connectivity
    void activate_overset();

    //! Activate mesh mapping
    void activate_mesh_map();

    bool has_mesh_mapping() const { return m_mesh_mapping; }

    void set_during_overset_advance(const bool flag)
    {
        m_during_overset_advance = flag;
    }

    bool is_during_overset_advance() const { return m_during_overset_advance; }

    bool is_anelastic() const;

private:
    amrex::AmrCore& m_mesh;

    SimTime m_time;

    mutable FieldRepo m_repo;

    pde::PDEMgr m_pde_mgr;

    PhysicsMgr m_physics_mgr;

    FieldBoundaryMgr m_field_boundary_mgr;

    std::unique_ptr<transport::TransportModel> m_transport;

    std::unique_ptr<turbulence::TurbulenceModel> m_turbulence;

    std::unique_ptr<IOManager> m_io_mgr;

    std::unique_ptr<PostProcessManager> m_post_mgr;

    std::unique_ptr<OversetManager> m_overset_mgr;

    std::unique_ptr<MeshMap> m_mesh_map;

    std::unique_ptr<ExtSolverMgr> m_ext_solver_mgr;

    std::unique_ptr<HelicsStorage> m_helics;

    bool m_mesh_mapping{false};

    // State of solver - know if during an overset timestep or not
    bool m_during_overset_advance{false};

    MultiBlockContainer* m_mbc{nullptr};
    ReadERFFunction m_read_erf{nullptr};
};

} // namespace kynema_sgf

#endif /* CFDSIM_H */
