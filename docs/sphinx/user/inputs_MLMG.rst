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
   include: smoother, bicgstab, cg, bicgcg, cgbicg, hypre, and petsc. 
   The hyper and petsc options will require compiling with those libraries.

.. input_param:: diffusion.bottom_verbose

   **type:** Integer, optional, default = 0

   Sets the verbosity of the bottom solver within MLMG.

.. input_param:: diffusion.bottom_rtol

   **type:** Real, optional, default = 1.0e-4

   Set the relative tolerance for the bottom solver for convergence.

.. input_param:: diffusion.bottom_atol

   **type:** Real, optional, default = -1.0

   Set the absolute tolerance for the bottom solve. Setting a negative number
   disables absolute tolerance check.

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



**MLMG performance tuning options**

The following options control how MLMG builds and coarsens the multigrid
hierarchy. They can be used to reduce the cost of the linear solves,
particularly at high MPI rank counts. As with the options above, the prefix
selects the solve; the same options apply to each. In the starting values
below, replace ``<solve>`` with each solve prefix you want to tune, i.e.
``diffusion``, ``nodal_proj``, and ``mac_proj`` (or a specific solve such as
``temperature_diffusion``, ``sdr``, or ``tke``)::

   <solve>.mg_rtol                  = 1.0e-6   # 1.0e-7 for the projections
   <solve>.mg_atol                  = 1.0e-11
   <solve>.deterministic            = false
   <solve>.do_agglomeration         = true
   <solve>.do_consolidation         = true
   <solve>.do_semicoarsening        = false    # set true only if grids are anisotropic
   <solve>.agg_grid_size            = 16        # try 32 on high rank counts
   <solve>.con_grid_size            = 16        # try 32 on high rank counts
   <solve>.max_coarsening_level     = 100       # coarsen as deep as possible
   <solve>.max_semicoarsening_level = 30
   <solve>.bottom_rtol              = 1.0e-4    # loosen if the bottom solve dominates
   <solve>.bottom_atol              = -1.0      # negative disables the absolute check

For example, to apply these to all three solves you would set
``diffusion.mg_rtol``, ``nodal_proj.mg_rtol``, and ``mac_proj.mg_rtol``, and so
on for each option.

The solver tolerances :input_param:`diffusion.mg_rtol` and
:input_param:`diffusion.mg_atol` are documented above under **MLMG options**.

The maximum coarsening level, :input_param:`diffusion.max_coarsening_level`,
is documented above under **Linear operator options**. Setting it to a large
value (e.g. ``100``) lets MLMG coarsen as deep as the grid allows, which is
generally desirable for performance.

The relative and absolute bottom-solver tolerances,
:input_param:`diffusion.bottom_rtol` and :input_param:`diffusion.bottom_atol`,
are documented above under **Bottom solver options**. Loosening
``bottom_rtol`` can reduce cost when the bottom solve dominates the total
solve time.

.. input_param:: diffusion.deterministic

   **type:** Boolean, optional, default = false

   If ``true``, the solver produces bit-for-bit reproducible results regardless
   of the number of ranks or the box distribution, at the cost of extra
   communication. Leave ``false`` for best performance; set ``true`` only when
   reproducibility is required (e.g. regression testing).

.. input_param:: diffusion.do_agglomeration

   **type:** Boolean, optional, default = true

   If ``true``, coarse-level grids are agglomerated (merged) onto fewer, larger
   boxes as the hierarchy is coarsened. This reduces the number of small boxes
   and the associated communication overhead on coarse levels.

.. input_param:: diffusion.do_consolidation

   **type:** Boolean, optional, default = true

   If ``true``, coarse levels are consolidated onto a subset of ranks so that
   fewer ranks participate in the coarsest solves. This lowers communication
   cost at high rank counts.

.. input_param:: diffusion.do_semicoarsening

   **type:** Boolean, optional, default = false

   If ``true``, MLMG may coarsen in only a subset of directions
   (semi-coarsening) rather than in all directions at once. Useful when the
   grid is anisotropic (very different cell sizes or extents per direction);
   leave ``false`` for roughly isotropic grids.

.. input_param:: diffusion.agg_grid_size

   **type:** Integer, optional, default = -1

   Target box size used when agglomerating coarse grids (see
   :input_param:`diffusion.do_agglomeration`). A value of ``16`` is a good
   starting point; ``32`` can perform better at high rank counts. A negative
   value lets AMReX choose.

.. input_param:: diffusion.con_grid_size

   **type:** Integer, optional, default = -1

   Target box size used when consolidating coarse grids (see
   :input_param:`diffusion.do_consolidation`). A value of ``16`` is a good
   starting point; ``32`` can perform better at high rank counts. A negative
   value lets AMReX choose.

.. input_param:: diffusion.max_semicoarsening_level

   **type:** Integer, optional, default = 0

   Maximum number of semi-coarsening levels allowed when
   :input_param:`diffusion.do_semicoarsening` is ``true``. A larger value (e.g.
   ``30``) allows the anisotropic coarsening to proceed as deep as possible.

.. input_param:: diffusion.semicoarsening_direction

   **type:** Integer, optional, default = -1

   Restricts semi-coarsening to a specific direction (``0``, ``1``, or ``2``)
   when :input_param:`diffusion.do_semicoarsening` is ``true``. Useful for
   strongly anisotropic grids where the fine direction is known. The default of
   ``-1`` lets AMReX select the direction automatically.
