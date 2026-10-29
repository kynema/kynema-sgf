.. _inputs_mlmg:

Section: MLMG options
~~~~~~~~~~~~~~~~~~~~~

.. tip::

   For further information the user is directed to the `AMReX linear
   solver documentation
   <https://amrex-codes.github.io/amrex/docs_html/LinearSolvers_Chapter.html>`_
   and other `AMReX runtime parameters
   <https://amrex-codes.github.io/amrex/docs_html/RuntimeParameters.html>`_


This section specifies the Multi-Level Multi-Grid (MLMG) options for each type
of linear solve. There are three types of linear solves performed in Kynema-SGF
"diffusion" which is a cell based Helmholtz like solve to advance the momentum
equations, "nodal_proj" is a node based pressure projection, and "mac_proj"
projects velocities to faces. The options are the same for each and the prefix
determines which MLMG option is being specified. Below the diffusion options are
described but the same options apply to "nodal_proj" and "mac_proj". It is also
possible to specify diffusion solver options for specific equations such as
temperature, to do that use "temperature_diffusion" as your prefix.


**Linear operator options**
   
.. input_param:: diffusion.max_coarsening_level

   **type:** Integer, optional, default = 100
   
   This parameter sets the max number of multigrid coarsening allowed at the lowest amr 
   level to the solver. 
   Typically setting to a large number means coarsen as much as possible until the grid 
   can not be coarsened anymore.
      
.. input_param:: diffusion.max_order

   **type:** Integer, optional, default = 2
   
   Order of the one-sided stencil applied near physical boundaries and fine/coarse boundaries.

**MLMG options**

.. input_param:: diffusion.multigrid_type

   **type:** String, optional, default = ``geometric``

   Sets how MLMG solves the coarsest AMR level (level 0). Valid choices are:

   - ``geometric``: geometric multigrid V-cycles (or F-cycles, see
     :input_param:`diffusion.fmg_maxiter`) on the coarsened levels.
   - ``algebraic``: the AMReX algebraic multigrid solver (AlgMG) solves level 0
     in every MLMG iteration. This can reduce the iteration count on meshes
     with a high cell aspect ratio, where geometric coarsening with a point
     smoother is less effective (see also ``diffusion.do_semicoarsening``).
     The coarsened geometric levels are not used for level 0, so setting
     :input_param:`diffusion.max_coarsening_level` to ``0`` avoids building
     them. The AlgMG solve of level 0 uses the bottom solver tolerances
   (:input_param:`diffusion.bottom_rtol` and
   :input_param:`diffusion.bottom_atol`). The bottom solver type is not
     used for level 0.
   - ``hybrid``: geometric until the iteration stalls or diverges, then
     algebraic, restarting from the best iterate so far. This requires
     :input_param:`diffusion.do_fixed_iters` ``= false``, because the switch
     needs a convergence test.

   Finer AMR levels always use geometric multigrid. ``algebraic`` and
   ``hybrid`` require a single-component operator. They can be used for
   ``mac_proj``, ``nodal_proj``, scalar diffusion solves (for example
   ``temperature_diffusion``) and the velocity diffusion solve with
   ``velocity_diffusion.use_tensor_operator = false`` and
   ``velocity_diffusion.use_segregated_operator = true``. The default velocity
   diffusion solve has three components. It reads ``diffusion`` and then
   ``velocity_diffusion``, so when ``diffusion.multigrid_type`` is not
   ``geometric``, also set ``velocity_diffusion.multigrid_type = geometric``.
   ``algebraic`` and ``hybrid`` also require ``diffusion.do_nsolve = false``
   (the default). In builds with ``KYNEMA_SGF_ENABLE_FFT``, the MAC
   projection uses FFT when ``mac_proj.use_fft = true`` (the default) and the
   case allows it; ``mac_proj.multigrid_type`` then has no effect.
   Invalid values and combinations abort with a message naming the input.

.. input_param:: diffusion.verbose

   **type:** Integer, optional, default = 0

   Sets the verbosity of the MLMG solver.

.. input_param:: diffusion.maxiter

   **type:** Integer, optional, default = 200

   Sets the max number of multigrid iterations. If :input_param:`diffusion.do_fixed_iters`
   is set to True, then AMReX will not abort if specified tolerance is not met
   after max iterations, otherwise it will abort.

.. input_param:: diffusion.do_fixed_iters

   **type:** Boolean, optional, default = true

   If ``true``, then AMReX will not abort if the specified tolerance is not met
   even after :input_param:`diffusion.maxiter` iterations have completed.

.. input_param:: diffusion.mg_rtol

   **type:** Real, optional, default = 1.0e-11
   
   Set the relative tolerance for the linear solver
   
.. input_param:: diffusion.mg_atol

   **type:** Real, optional, default = 1.0e-14
   
   Set the absolute tolerance for the linear solver

.. input_param:: diffusion.fmg_maxiter

   **type:** Integer, optional, default = 0

   Sets the number of F-cycle MG iterations to perform before switching to V-cycle MG.

.. input_param:: diffusion.num_pre_smooth

   **type:** Integer, optional, default = 2

   Number of pre smoothing steps

.. input_param:: diffusion.num_post_smooth

   **type:** Integer, optional, default = 2

   Number of post smoothing steps

.. input_param:: diffusion.num_final_smooth

   **type:** Integer, optional, default = 8

   Number of final smoother steps applied

.. input_param:: diffusion.num_bottom_smooth

   **type:** Integer, optional, default = 0

   Number of smoother steps applied during bottom solve.

**Bottom solver options**
   
.. input_param:: diffusion.bottom_solver

   **type:** String, optional, default = "bicgstab"
   
   Set the bottom solver type. Current bottom solver options 
   include: smoother, bicgstab, cg, bicgcg, cgbicg, amg, hypre, and petsc.
   The amg option uses the AMReX algebraic multigrid solver (AlgMG) as the
   bottom solver.
   The hyper and petsc options will require compiling with those libraries.

.. input_param:: diffusion.bottom_verbose

   **type:** Integer, optional, default = 0

   Sets the verbosity of the bottom solver within MLMG.

.. input_param:: diffusion.bottom_rtol

   **type:** Real, optional, default = 2.2e-4 (double precision), 1.2e-4
   (single precision)

   Set the relative tolerance for the bottom solver for convergence.

.. input_param:: diffusion.bottom_atol

   **type:** Real, optional, default = 2.2e-12 (double precision), 1.2e-5
   (single precision); -1.0 when :input_param:`diffusion.multigrid_type` is
   ``algebraic`` or ``hybrid``

   Set the absolute tolerance for the bottom solve. Setting a negative number
   disables absolute tolerance check. With ``algebraic`` or ``hybrid``, the
   bottom tolerances also control the AlgMG solve of level 0. A positive
   absolute tolerance can then stop that solve before the MLMG tolerances are
   met, and the MLMG iterations stall. So the check is disabled unless this
   input is set.

.. input_param:: diffusion.bottom_maxiter

   **type:** Integer, optional, default = 200

   Maximum number of iterations for the bottom solver

.. input_param:: diffusion.hypre_interface

   **type:** String, optional, default = ``ij``

   The hypre interface to use when :input_param:`diffusion.bottom_solver` is set
   to ``hypre``. Valid choices are: ``ij``, ``semi_structured``, and
   ``structured``.

.. input_param:: diffusion.hypre_namespace

   The ParmParse ``prefix`` where the hypre options must be read from for this
   solver. For example, to set hypre options for NodalProjector

   ..
      nodal_proj.hypre_namespace = "nodal_proj.hypre"
      nodal_proj.hypre.hypre_solver = GMRES
      nodal_proj.hypre.hypre_preconditioner = BoomerAMG


