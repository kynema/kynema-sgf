.. _inputs_amr:

Section: AMReX and AMR
~~~~~~~~~~~~~~~~~~~~~~

.. tip::

   The user is encouraged to refer to the full list of `AMReX runtime
   parameters
   <https://amrex-codes.github.io/amrex/docs_html/RuntimeParameters.html>`_
   for more advanced usage. Another section of particular interest is
   the `description of AMReX grid creation
   <https://amrex-codes.github.io/amrex/docs_html/GridCreation.html#sec-grid-creation>`_.

This section contains some of the common input parameters used by the
core AMReX mesh data structure ``AmrCore`` to determine the base mesh
and adaptive mesh refinement strategies. The refinement criteria is
specified through :ref:`inputs_tagging`.

.. input_param:: amr.n_cell

   **type:** List of 3 integers, mandatory

   The number of cells in the coarsest level of AMR hierarchy

.. input_param:: amr.max_level

   **type:** Integer, mandatory

   The maximum AMR level in the refinement hierarchy.

.. input_param:: amr.max_grid_size

   **type:** Integer, optional, default: [build dependent]

   Maximum number of cells at level 0 in each grid.
   There are also options to specify this value in each direction,
   please refer to AMReX documentation.

.. input_param:: amr.blocking_factor

   **type:** Integer, optional, default: 8

   Each grid must be divisible by :input_param:`amr.blocking_factor`.
   There are also options to specify this value in each direction,
   please refer to AMReX documentation.



