.. _inputs_immersedterrain:

Section: ImmersedTerrain
~~~~~~~~~~~~~~~~~~~~~~~~

These parameters are active when ``ImmersedTerrain`` is included in
:input_param:`incflo.physics`. ImmersedTerrain is the successor to
:ref:`TerrainDrag <inputs_terraindrag>`: instead of a binary blanking it
stores the fraction of each cell occupied by terrain, so slopes are not
represented as a staircase. It is used together with the
``ImmersedDragForcing`` momentum source (see
:ref:`inputs_momentum_sources`).

ImmersedTerrain declares the following fields, each with one ghost cell:

- ``terrain_fraction`` (Real): fraction of the cell volume occupied by terrain,
  0 in fluid, 1 inside terrain, in between for cells cut by the surface: the
  fraction of the cell column below the terrain height at the cell center.
- ``terrain_mask`` (Int): 0 for fluid, 1 for cells fully inside the terrain,
  2 for surface cells, i.e. partially filled cells and fluid cells that share
  a face with a solid cell on any of the six sides. The mask can be used
  directly with
  ``FieldRefinement``: ``field_error = 1.5`` tags only the surface band and
  ``field_error = 0.5`` tags surface and solid cells.
- ``terrain_surface`` (Real): terrain height at the cell center.
- ``terrain_diffusion_xf/yf/zf`` (Real, face centered): factors applied to the
  face diffusion coefficients so that the flux across a wall face is the
  no-slip flux :math:`\mu u / d_1` at the true wall distance :math:`d_1`. With
  ``center`` weighting :math:`d_1 = z_k - h` for the first cell whose center is
  above the terrain; with ``fraction`` weighting the partial cell is the wall
  cell and :math:`d_1 = (1-\beta)\Delta z/2`, the height of the centroid of
  its fluid part, and the face above it is scaled by :math:`2/(2-\beta)`.
  Side walls use :math:`d_1 = \Delta_f/2`. The factors apply to every
  equation solved with the shared diffusion operator.

.. input_param:: ImmersedTerrain.terrain_file

   **type:** String, optional, default = ``terrain.amrwind``

   Input file for terrain height data, in the same flat-grid format as
   :input_param:`TerrainDrag.terrain_file`.

.. input_param:: ImmersedTerrain.implicit_projection

   **type:** Boolean, optional, default = false

   Apply the immersed drag implicitly through the nodal and MAC projections
   instead of as an explicit source term. The terrain then behaves as a fluid of
   density :math:`\rho (1 + C \Delta t)` in the pressure solve, with
   :math:`C = w C_d / \Delta z`, the drag weight :math:`w` of
   :input_param:`ImmersedTerrain.drag_weight` and :math:`C_d` taken from
   :input_param:`ImmersedDragForcing.drag_coefficient`, so that the pressure
   gradient produces no velocity inside the terrain. Without it the projection
   re-injects :math:`\Delta t \nabla p / \rho` inside the body every step and the
   velocity residual inside the terrain decreases only linearly with the time
   step. When active, ``ImmersedDragForcing`` adds no explicit drag. Declares
   the field ``terrain_drag_rate``.

.. input_param:: ImmersedTerrain.drag_weight

   **type:** String, optional, default = ``fraction``

   Weight of the immersed drag in partially filled cells. ``fraction`` uses the terrain
   fraction :math:`\beta`. ``center`` treats a cell whose center is inside the terrain
   (:math:`\beta \ge 0.5`) as fully solid, with a wall on the faces next to it, and
   any other partial cell as a fluid cell with no drag. With laminar flow and ``fraction``,
   the partial cell takes no drag either: it is a fluid
   cell whose wall is carried by the no-slip flux at its fluid centroid, and only cells
   entirely inside the terrain are held by the drag. Applies to the explicit drag and to
   the implicit projection rate.
