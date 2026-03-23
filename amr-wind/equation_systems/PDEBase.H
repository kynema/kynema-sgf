#ifndef PDEBASE_H
#define PDEBASE_H

#include <string>

#include "amr-wind/core/Factory.H"
#include "amr-wind/core/CollMgr.H"
#include "amr-wind/equation_systems/PDEFields.H"
#include "amr-wind/core/FieldDescTypes.H"
#include "amr-wind/incflo_enums.H"

namespace amr_wind {

class CFDSim;

namespace pde {

/**
 *  \defgroup eqsys Equation Systems
 *
 *  Governing partial-differential equations
 *
 *  The governing equations and their solution procedure is implemented within
 *  equation systems that all inherit from PDEBase class. Each term of the
 *  governing equation is implemented as a \ref pdeop that can be specialized
 *  for different equations. The behavior is specialized using two traits: 1. a
 *  PDE trait, e.g., ICNS representing the incompressible Navier-Stokes
 *  equation, or the Density trait representing the scalar density transport
 *  equation; and 2. a Scheme trait, e.g., Godunov or MOL that determines the
 *  discretization process.
 */

/** Abstract interface defining a PDE and its solution process withn AMR-Wind
 *  \ingroup eqsys
 *
 *  The exact behavior is implemented by subclass PDESystem that is specialized
 *  using two traits: the PDE type trait, and the numerical scheme trait.
 *  Currently, two numerical scheme traits are supported: Godunov and MOL.
 *
 *  Each PDESystem instance contains its own PDEFields object that holds fields
 *  unique to the PDE. These fields are: the variable solved by the PDE (e.g.,
 *  velocity, temperature, etc.), the effective dynamic viscosity, the source
 *  terms, the diffusion term, and the convective term. A PDE may declare
 *  additional fields that are not stored in this object, but can be queried
 *  using the FieldRepo instance. For example, the incompressible Navier-Stokes
 *  (ICNS) always declares density and pressure in addition to the velocity
 *  fields.
 *
 */
class PDEBase : public Factory<PDEBase, CFDSim&>
{
public:
    ~PDEBase() override = default;

    //! Return the object that holds references to fields for this PDE system
    virtual PDEFields& fields() = 0;
    [[nodiscard]] virtual const PDEFields& fields() const = 0;

    //! Perform initialization actions for a PDE after the mesh is generated
    virtual void initialize() = 0;

    //! Perform updates specific to PDE after a regrid is performed
    virtual void post_regrid_actions() = 0;

    //! Compute the source (forcing) terms for this PDE
    virtual void compute_source_term(FieldState fstate) = 0;

    //! Compute the effective dynamic viscosity for this PDE
    virtual void compute_mueff(FieldState fstate) = 0;

    //! Compute the diffusion term used in the RHS of the PDE system
    virtual void compute_diffusion_term(FieldState fstate) = 0;

    //! Perform necessary steps to prepare for advection calculations
    virtual void pre_advection_actions(FieldState fstate) = 0;

    //! Compute the time derivative and advective term for the PDE system
    virtual void compute_advection_term(FieldState fstate) = 0;

    /** Combine the source, diffusion, and advection term to obtain RHS
     *
     *  This method behaves differently depending upon whether the diffusion
     *  term is treated implicitly, with Crank-Nicolson scheme, or explicitly.
     *  The predictor step only uses the "old" states for RHS computation.
     */
    virtual void compute_predictor_rhs(DiffusionType difftype) = 0;

    /** Combine the source, diffusion, and advection term to obtain RHS
     *
     *  This method behaves differently depending upon whether the diffusion
     *  term is treated implicitly, with Crank-Nicolson scheme, or explicitly.
     *  The corrector step uses both "old" state as well as the updated "new"
     *  state from predictor step to compute the RHS term.
     */
    virtual void compute_corrector_rhs(DiffusionType difftype) = 0;

    //! Solve the diffusion linear system and update the field
    virtual void solve(amrex::Real dt) = 0;

    //! Perform post-processing actions after a system solve
    virtual void post_solve_actions() = 0;

    virtual void improve_explicit_diffusion(amrex::Real dt) = 0;

    //! Base class identifier used for factory registration interface
    static std::string base_identifier() { return "PDESystem"; }
};

/** PDE systems manager
 *
 *  This class manages a collection of PDEs and provides an interface to
 *  register PDEs. An incompressible N-S equation is always registered and
 *  treated as a special case. The remaining PDEs are created on-demand and are
 *  accessed through the equation systems vector.
 */
class PDEMgr : public CollMgr<PDEMgr, PDEBase>
{
public:
    explicit PDEMgr(CFDSim& sim);

    ~PDEMgr() = default;

    //! Return the incompressible Navier-Stokes instance
    PDEBase& icns() { return *m_icns; }
    [[nodiscard]] const PDEBase& icns() const { return *m_icns; }

    //! Create the incompressible Navier-Stokes instance and return its
    //! reference
    PDEBase& register_icns();

    //! Register a new PDE instance and return its reference
    PDEBase& register_transport_pde(const std::string& pde_name);

    //! Advance states for all registered PDEs at the beginning of a timestep
    void advance_states();

    //! Prepare boundaries for registered PDEs at the beginning of a timestep:
    //! -- intended for when boundaries rely on intermediate time states (nph)
    void prepare_boundaries();

    //! Check that the bounds of density correspond to the expected values
    //! -- only applies to two-phase and constant-density single-phase
    //! -- performed once in post_init, will abort if density is unexpected
    void density_check();

    //! Call fillpatch operator on state variables for all registered PDEs
    void fillpatch_state_fields(
        amrex::Real time, FieldState fstate = FieldState::New);

    /** Return the vector containing all registered PDE instances
     *
     *  Note that this does not contain the ICNS system as it is treated
     *  separately for algorithmic reasons.
     */
    TypeVector& scalar_eqns() { return m_obj_vec; }
    [[nodiscard]] const TypeVector& scalar_eqns() const { return m_obj_vec; }

    //! Query if a PDE has already been registered
    [[nodiscard]] bool has_pde(const std::string& pde_name) const;

    [[nodiscard]] int num_ghost_state() const;

    [[nodiscard]] bool constant_density() const { return m_constant_density; }

    [[nodiscard]] std::string scheme() const { return m_scheme; }

private:
    //! Instance of the CFD simulation controller
    CFDSim& m_sim;

    //! ICNS instance
    std::unique_ptr<PDEBase> m_icns;

    //! Unique identifier for the numerical scheme
    std::string m_scheme;

    //! Flag indicating whether Godunov scheme is active
    bool m_use_godunov{true};

    //! Flag indicating whether density is constant for this simulation
    bool m_constant_density{true};
};

} // namespace pde
} // namespace amr_wind

#endif /* PDEBASE_H */
