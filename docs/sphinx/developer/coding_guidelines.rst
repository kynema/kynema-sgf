.. dev_coding-guidelines:

Coding Guidelines
=================

This section presents an incomplete list of coding guidelines currently used to
develop Kynema-SGF code. As is the case with rules, it might not apply to all
circumstances. Please use your best judgment when trying to comply with the
rules and break them if necessary, but be ready to convince the rest of the
team.

**Use modern C++**

   Traditionally, AMReX based solvers have used C++ for managing data structures
   but Fortran for computational kernels. However, Kynema-SGF is purely written in
   C++. One of the primary drivers for this choice is our desire to target
   next-generation Exascale systems that require targeting GPUs through CUDA,
   HIP, and other offloading paradigms. Support for languages other than C/C++
   is not very mature for these architectures.

   Consult `ISO C++ Core Guidelines
   <https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines>`_ for best
   practices when coding in C++.

**Write code targeting heterogeneous architectures**

   All code must be written in a way such that it can run on GPUs without
   errors. Please consult the `GPU section
   <https://amrex-codes.github.io/amrex/docs_html/GPU_Chapter.html>`_ in AMReX
   documentation to learn the correct data structures and looping paradigms to
   use when targeting heterogeneous architectures.

**Unit and regression tests**

   Code modifications must not break existing unit tests or regression tests. If
   a change is known to break current behavior, e.g., a bug fix or an
   enhancement, the tests must be updated to account for this change and any
   necessary updates to documentation must be made to document this.

   New features must be accompanied with the necessary suite of unit and
   regression tests as appropriate. There should be at least one test that
   exercises the feature to catch issue early on. New features must also be
   documented in the user and theory manual.

**Pull-request based workflow**

   All updates must be submitted as `pull-requests
   <https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/proposing-changes-to-your-work-with-pull-requests>`_
   on the public GitHub repository. Pull requests will be reviewed by core
   developers and must also pass all Continuous Integration checks before it is
   merged into mainline source tree.

   If a pull-requests fails during nightly tests after it is merged with the
   mainline code base, it will be reverted and a new pull request should be
   submitted with appropriate fixes.


Style Guide/Recommendations
----------------------------

This section documents the preferred style for formatting source code within
Kynema-SGF. While the guidelines presented in the previous section are based on
technical and quality assurance reasons, the formatting style is largely a
matter of individual taste and there is no hard and fast rule. However, to
ensure consistency in a codebase used by a large team, we provide some
guidelines here.

Kynema-SGF comes with a :file:`.clang-format` definition that can be used with
`Clang format <https://clang.llvm.org/docs/ClangFormat.html>`_ to automatically
reformat source code to be consistent with the rest of the codebase. Many source
code editors come with Clang format integration that allows you to reformat code
from within your text editor. Please consult the Clang format documentation, or
your editor's documentation to setup clang-format integration. However, please
be careful that you only use Clang format on code you've written modified and
not on code that has not been written by you.

Other Kynema-SGF specific conventions:

#. All Kynema-SGF specific code must be within ``kynema_sgf`` namespace. Code for
   unit tests must be in ``kynema_sgf_tests`` namespace.

#. Following AMReX convention, header files will use a ``.H`` extension and C++
   source files will use a ``.cpp`` extension.

#. Use `snake_case <https://en.wikipedia.org/wiki/Snake_case>`_ almost
   everywhere, i.e., for variables, namespaces, function and method names.
   Exceptions include when overriding AMReX class methods in inherited classes.

#. Use `CamelCase <https://en.wikipedia.org/wiki/Camel_case>`_ for class names.
   Capitalize the first letter of the first word in the compound name also.

#. Use ``m_`` prefix for class instance variables. No prefix for class methods.

#. Keep logic in functions short. Always use descriptive names for function
   names and variables that provide the reader with clues as to what the
   function does.

Sample C++ code
~~~~~~~~~~~~~~~

.. code-block:: c++
   :linenos:

   #ifndef CFDSIM_H
   #define CFDSIM_H

   #include "AMReX_AmrCore.H"
   #include "SimTime.H"
   #include "FieldRepo.H"
   #include "PDEBase.H"
   #include "Physics.H"

   namespace kynema_sgf {

   // Forward declaration if necessary
   namespace turbulence {
   class TurbulenceModel;
   }

   /** Data structures for a CFD simulation
    *
    *  CFDSim is a thin wrapper that holds all the necessary objects for a CFD
    *  simulation. The key data members within this object are:
    *
    *  - mesh (amrex::AmrCore) The AMR mesh hierarchy data structure
    *  - time (SimTime)        The time object
    *  - repo (FieldRepo)      The field repository
    *  - pde_manager (PDEMgr)  PDE manager interface
    *  - physics_manager       List of physics active in this simulation
    *  - turbulence_model      Reference to the turbulence model
    */
   class CFDSim
   {
   public:
       CFDSim(amrex::AmrCore& mesh);

       ~CFDSim();

       //! Return the AMR mesh hierarchy
       amrex::AmrCore& mesh() { return m_mesh; }
       const amrex::AmrCore& mesh() const { return m_mesh; }

       //! Return simulation time control
       SimTime& time() { return m_time; }
       const SimTime& time() const { return m_time; }

       //! Return the field repository
       FieldRepo& repo() { return m_repo; }
       const FieldRepo& repo() const { return m_repo; }

       //! Return the PDE manager instance
       pde::PDEMgr& pde_manager() { return m_pde_mgr; }
       const pde::PDEMgr& pde_manager() const { return m_pde_mgr; }

       //! Return the physics manager instance
       PhysicsMgr& physics_manager() { return m_physics_mgr; }
       const PhysicsMgr& physics_manager() const { return m_physics_mgr; }

       //! Return a vector of physics instances active in this simulation
       PhysicsMgr::TypeVector& physics() { return m_physics_mgr.objects(); }
       const PhysicsMgr::TypeVector& physics() const { return m_physics_mgr.objects(); }

       //! Return the turbulence model instance used in this simulation
       turbulence::TurbulenceModel& turbulence_model() { return *m_turbulence; }
       const turbulence::TurbulenceModel& turbulence_model() const
       { return *m_turbulence; }

       //! Initialize the turbulence model after reading necessary inputs
       void create_turbulence_model();

       //! Initialize the different physics models based on user inputs
       void init_physics();

   private:
       amrex::AmrCore& m_mesh;

       SimTime m_time;

       FieldRepo m_repo;

       pde::PDEMgr m_pde_mgr;

       PhysicsMgr m_physics_mgr;

       std::unique_ptr<turbulence::TurbulenceModel> m_turbulence;
   };

   } // namespace kynema_sgf

   #endif /* CFDSIM_H */
