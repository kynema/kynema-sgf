#ifndef FIELDREPO_H
#define FIELDREPO_H

#include <string>
#include <unordered_map>

#include "amr-wind/core/FieldDescTypes.H"
#include "amr-wind/core/FieldUtils.H"
#include "amr-wind/core/Field.H"
#include "amr-wind/core/IntField.H"
#include "amr-wind/core/ScratchField.H"
#include "amr-wind/core/IntScratchField.H"

#include "AMReX_AmrCore.H"
#include "AMReX_MultiFab.H"
#include "AMReX_iMultiFab.H"

namespace amr_wind {

/**
 *  \defgroup fields Field management
 *  Field management infrastructure
 *
 *  Field is at the core of AMR-Wind implementation that handles all the data
 *  used during the simulation. FieldRepo manages a collection of fields and is
 *  responsible for all memory allocation and reinitialization during regrid
 *  process. Field provides a higher-level interface than MultiFabs and also
 *  provides a unified interface to manage boundary conditions and states.
 *
 *  \ingroup core
 */

/** Object that holds the MultiFab instances at all levels for a given field
 */
struct LevelDataHolder
{
    LevelDataHolder();

    //! real multifabs for all the known fields at this level
    amrex::Vector<amrex::MultiFab> m_mfabs;
    //! Factory for creating new FABs
    std::unique_ptr<amrex::FabFactory<amrex::FArrayBox>> m_factory;

    //! int fabs for all known fields at this level
    amrex::Vector<amrex::iMultiFab> m_int_fabs;
    std::unique_ptr<amrex::FabFactory<amrex::IArrayBox>> m_int_fact;
};

/** Field Repository
 *  \ingroup fields
 *
 *  FieldRepo is a database that manages handling of fields used during a
 *  simulation. It provides an interface to automatically handle the
 *  creation/destruction of lower-level data structures (amrex::MultiFab) during
 *  regrids and handles the interface between the actual data and the
 *  amr_wind::Field interface which doesn't hold any actual data.
 *
 *  FieldRepo allows the user to declare instances of fields and their types
 *  (e.g., cell or node cenrered), and allows querying and accessing these
 *  fields at different places in the code and ensures the same data is accessed
 *  everywhere.
 *
 *  Fields can be declared at any time, even before the actual mesh instance is
 *  created, and it is recommended that they be created before mesh is
 *  initialized so that initialization happens automatically. When fields are
 *  declared before a mesh is available, data allocation is deferred until the
 *  mesh is generated (`MakeNewLevelFromScratch`). However, if a mesh is
 *  available, data is allocated and ready for use immediately upon declaration.
 *
 *  Field declarations use a unique string name to differentiate different
 *  fields. If amr_wind::FieldRepo::declare_field is called multiple times with
 *  the same field name, the previously created instance is returned on
 *  subsequent calls. However, an error is raised if the subsequent calls do not
 *  have the same field parameters (e.g., number of components, states, or ghost
 *  cells).
 *
 *  Fields can be accessed using amr_wind::FieldRepo::get_field method. Unlike
 *  the declaration method, this method will raise an error if no field exists.
 *  amr_wind::FieldRepo::field_exists can be used to determine if a field exists
 *  in the repository.
 *
 *  FieldRepo also manages integer fields (IntField) as well as creation of
 *  ScratchField instances.
 */
class FieldRepo
{
public:
    friend class Field;
    friend class IntField;

    explicit FieldRepo(const amrex::AmrCore& mesh)
        : m_mesh(mesh), m_leveldata(mesh.maxLevel() + 1)
    {}

    FieldRepo(const FieldRepo&) = delete;
    FieldRepo& operator=(const FieldRepo&) = delete;
    ~FieldRepo() = default;

    //! Perform field data management tasks when
    //! amrex::AmrCore::MakeNewLevelFromScratch is called
    void make_new_level_from_scratch(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm);

    /** Allocate and initialize data on a finer level from the next coarser
     * level
     *
     *  This is called by the AmrCore automatically during regrid process
     */
    void make_new_level_from_coarse(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm);

    /** Allocate and initialize data on an existing level during regrid
     *
     *  This is called by the AmrCore automatically during regrid process
     */
    void remake_level(
        int lev,
        amrex::Real time,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm);

    //! Remove a level during regrid
    void clear_level(int lev);

    /** Create a new field instance
     *
     *  \param name Unique identifier for this field
     *  \param ncomp Number of components in this field (default: 1)
     *  \param ngrow Number of ghost cells/nodes for this field (default: 1)
     *  \param nstates Number of time states for this field (default: 1)
     *  \param floc Field location (default: cell-centered)
     */
    Field& declare_field(
        const std::string& name,
        int ncomp = 1,
        int ngrow = 0,
        int nstates = 1,
        FieldLoc floc = FieldLoc::CELL);

    /** Declare a cell-centered field
     *
     *  This is a convenience wrapper around `declare_field`
     *
     *  \sa amr_wind::FieldRepo::declare_field
     */
    Field& declare_cc_field(
        const std::string& name,
        const int ncomp = 1,
        const int ngrow = 0,
        const int nstates = 1)
    {
        return declare_field(name, ncomp, ngrow, nstates, FieldLoc::CELL);
    }

    /** Declare a node-centered field
     *
     *  This is a convenience wrapper around `declare_field`
     *
     *  \sa amr_wind::FieldRepo::declare_field
     */
    Field& declare_nd_field(
        const std::string& name,
        const int ncomp = 1,
        const int ngrow = 0,
        const int nstates = 1)
    {
        return declare_field(name, ncomp, ngrow, nstates, FieldLoc::NODE);
    }

    /** Declare a face field on the faces with normal in x-direction
     *
     *  This is a convenience wrapper around `declare_field`
     *
     *  \sa amr_wind::FieldRepo::declare_field
     */
    Field& declare_xf_field(
        const std::string& name,
        const int ncomp = 1,
        const int ngrow = 0,
        const int nstates = 1)
    {
        return declare_field(name, ncomp, ngrow, nstates, FieldLoc::XFACE);
    }

    /** Declare a face field on the faces with normal in y-direction
     *
     *  This is a convenience wrapper around `declare_field`
     *
     *  \sa amr_wind::FieldRepo::declare_field
     */
    Field& declare_yf_field(
        const std::string& name,
        const int ncomp = 1,
        const int ngrow = 0,
        const int nstates = 1)
    {
        return declare_field(name, ncomp, ngrow, nstates, FieldLoc::YFACE);
    }

    /** Declare a face field on the faces with normal in z-direction
     *
     *  This is a convenience wrapper around `declare_field`
     *
     *  \sa amr_wind::FieldRepo::declare_field
     */
    Field& declare_zf_field(
        const std::string& name,
        const int ncomp = 1,
        const int ngrow = 0,
        const int nstates = 1)
    {
        return declare_field(name, ncomp, ngrow, nstates, FieldLoc::ZFACE);
    }

    /** Declare a vector field with components in the x,y,z faces
     *
     *  This is a convenience function around `declare_field`
     */
    amrex::Vector<Field*> declare_face_normal_field(
        const amrex::Vector<std::string>& names,
        const int ncomp = 1,
        const int ngrow = 0,
        const int nstates = 1)
    {
        AMREX_ASSERT(names.size() == AMREX_SPACEDIM);
        return {
            &declare_field(names[0], ncomp, ngrow, nstates, FieldLoc::XFACE),
            &declare_field(names[1], ncomp, ngrow, nstates, FieldLoc::YFACE),
            &declare_field(names[2], ncomp, ngrow, nstates, FieldLoc::ZFACE),
        };
    }

    /** Return a previously created field identified by name and time state
     */
    Field& get_field(
        const std::string& name, FieldState fstate = FieldState::New) const;

    /** Return a previously created field using a unique identifier
     */
    Field& get_field(const int field_id) const
    {
        return *m_field_vec.at(field_id);
    }

    /** Return a previously created mesh mapping field using field location
     */
    Field& get_mesh_mapping_field(FieldLoc floc) const;

    /** Return a previously created mesh mapping detJ using field location
     */
    Field& get_mesh_mapping_det_j(FieldLoc floc) const;

    //! Query if field uniquely identified by name and time state exists in
    //! repository
    bool field_exists(
        const std::string& name, FieldState fstate = FieldState::New) const;

    /** Declare an integer field
     *
     *  Integer fields are used for masks and iblanking with overset meshes
     *
     *  \param name [in] Unique identifier for the field
     *  \param ncomp Number of components in this field (default: 1)
     *  \param ngrow Number of ghost cells/nodes for this field (default: 1)
     *  \param nstates Number of time states for this field (default: 1)
     *  \param floc Field location (default: cell-centered)
     */
    IntField& declare_int_field(
        const std::string& name,
        int ncomp = 1,
        int ngrow = 0,
        int nstates = 1,
        FieldLoc floc = FieldLoc::CELL);

    //! Return a reference to an integer field
    IntField& get_int_field(
        const std::string& name, FieldState fstate = FieldState::New) const;

    /** Return a reference to a previously created integer field using unique
     * identifier
     */
    IntField& get_int_field(const int field_id) const
    {
        return *m_int_field_vec.at(field_id);
    }

    //! Query if an integer field exists
    bool int_field_exists(
        const std::string& name, FieldState fstate = FieldState::New) const;

    /** Create a scratch field
     *
     *  ScratchField is a temporary field used to compute and store intermediate
     *  quantities that are needed within a timestep. Scratch fields provide an
     *  API similar to Field, but do not contain multiple states. The underlying
     *  MultiFab data is not managed by FieldRepo and, therefore, ScratchField
     *  do not survive a regrid. This method returns a unique_ptr instance that
     *  is only valid within a timestep. It is not safe to hold a reference to
     *  the ScratchField object across timesteps.
     */
    std::unique_ptr<ScratchField> create_scratch_field(
        const std::string& name,
        int ncomp = 1,
        int nghost = 0,
        FieldLoc floc = FieldLoc::CELL) const;

    /** Create a scratch field
     *
     *  ScratchField is a temporary field used to compute and store intermediate
     *  quantities that are needed within a timestep. Scratch fields provide an
     *  API similar to Field, but do not contain multiple states. The underlying
     *  MultiFab data is not managed by FieldRepo and, therefore, ScratchField
     *  do not survive a regrid. This method returns a unique_ptr instance that
     *  is only valid within a timestep. It is not safe to hold a reference to
     *  the ScratchField object across timesteps.
     */
    std::unique_ptr<ScratchField> create_scratch_field(
        int ncomp = 1, int nghost = 0, FieldLoc floc = FieldLoc::CELL) const;

    std::unique_ptr<ScratchField> create_scratch_field_on_host(
        const std::string& name,
        int ncomp = 1,
        int nghost = 0,
        FieldLoc floc = FieldLoc::CELL) const;

    std::unique_ptr<ScratchField> create_scratch_field_on_host(
        int ncomp = 1, int nghost = 0, FieldLoc floc = FieldLoc::CELL) const;

    std::unique_ptr<IntScratchField> create_int_scratch_field_on_host(
        const std::string& name,
        int ncomp = 1,
        int nghost = 0,
        FieldLoc floc = FieldLoc::CELL) const;

    std::unique_ptr<IntScratchField> create_int_scratch_field_on_host(
        int ncomp = 1, int nghost = 0, FieldLoc floc = FieldLoc::CELL) const;

    //! Advance all fields with more than one timestate to the new timestep
    void advance_states();

    //! Return a reference to the underlying AMR mesh instance
    const amrex::AmrCore& mesh() const { return m_mesh; }

    //! Total number of levels currently active in the AMR mesh
    int num_active_levels() const { return m_mesh.finestLevel() + 1; }

    //! Number of fields registered in the database
    int num_fields() const { return static_cast<int>(m_field_vec.size()); }

    //! Return list of fields registered (for unit-test purposes only)
    const amrex::Vector<std::unique_ptr<Field>>& fields() const
    {
        return m_field_vec;
    }

    //! Return factory instance at a given level
    const amrex::FabFactory<amrex::FArrayBox>& factory(int lev) const
    {
        return *m_leveldata[lev]->m_factory;
    }

protected:
    /** Return the amrex::MultiFab instance for a field at a given level
     *
     *  \param fid Unique integer field identifier for this field
     *  \param lev AMR level
     */
    amrex::MultiFab& get_multifab(const unsigned fid, const int lev)
    {
        BL_ASSERT(lev <= m_mesh.finestLevel());
        return m_leveldata[lev]->m_mfabs[fid];
    }

    /** Return the integer fab instance for a field at a given level
     *
     *  \param fid Unique integer field identifier for this field
     *  \param lev AMR level
     */
    amrex::iMultiFab& get_int_fab(const unsigned fid, const int lev)
    {
        BL_ASSERT(lev <= m_mesh.finestLevel());
        return m_leveldata[lev]->m_int_fabs[fid];
    }

    //! Create a new state for a field
    Field& create_state(Field& field, FieldState fstate);

    //! Allocate field data for a single level outside of regrid
    void allocate_field_data(
        int lev,
        const Field& field,
        LevelDataHolder& level_data,
        const amrex::FabFactory<amrex::FArrayBox>& factory);

    //! Allocate field data at all levels
    void allocate_field_data(Field& field);

    //! Allocate data at a level during regrid
    void allocate_field_data(
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm,
        LevelDataHolder& level_data,
        const amrex::FabFactory<amrex::FArrayBox>& factory);

    void allocate_field_data(
        int lev, const IntField& field, LevelDataHolder& level_data);

    void allocate_field_data(const IntField& field);

    void allocate_field_data(
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm,
        LevelDataHolder& level_data,
        const amrex::FabFactory<amrex::IArrayBox>& factory);

    //! Reference to the mesh instance
    const amrex::AmrCore& m_mesh;

    //! Array (size: nlevels) of data holder that contains another array of
    //! MultiFabs for all fields at that level.
    amrex::Vector<std::unique_ptr<LevelDataHolder>> m_leveldata;

    //! References to field instances identified by unique integer
    mutable amrex::Vector<std::unique_ptr<Field>> m_field_vec;

    //! Reference to integer field instances identified by unique integer
    mutable amrex::Vector<std::unique_ptr<IntField>> m_int_field_vec;

    //! Map of field name to unique integer ID for lookups
    std::unordered_map<std::string, size_t> m_fid_map;

    //! Map of integer field name to unique integer ID for lookups
    std::unordered_map<std::string, size_t> m_int_fid_map;

    //! Flag indicating if mesh is available to allocate field data
    bool m_is_initialized{false};
};

} // namespace amr_wind

#endif /* FIELDREPO_H */
