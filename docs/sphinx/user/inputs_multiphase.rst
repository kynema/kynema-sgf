.. _inputs_multiphase:

Section: MultiPhase
~~~~~~~~~~~~~~~~~~~~

This section is for the MultiPhase physics module, which is required for flows that solve air and water together using the volume-of-fluid method.

.. input_param:: MultiPhase.density_fluid1

   **type:** Real, optional, default = 10.

   Density of the liquid phase, i.e. the one that the vof variable refers to. Recommended value for water is near 1000.

.. input_param:: MultiPhase.density_fluid2

   **type:** Real, optional, default = 1.

   Density of the gas phase, i.e. the one that the vof variable does not refer to. Recommended value for air is near 1.

.. input_param:: MultiPhase.verbose

   **type:** Integer, optional, default = 0

   Verbosity of the ``MultiPhase`` physics output. Verbosity greater than 0 will output conservation statistics, recording the
   total water volume and air volume in the domain, the differences between these and the initial volumes, and the difference
   between the initial momentum and the current momentum. These quantities can be used to confirm conservation properties
   in periodic cases without source terms.

.. input_param:: MultiPhase.water_level

   **type:** Real, optional, default = 0.

   Reference water level of the simulated flow. This value is used by the ``MultiPhase`` physics module to create
   a reference density field, which enables the use of a perturbational gravity term through the ICNS source
   term ``GravityForcing``. This parameter is also used to adjust how the ``GeostrophicForcing`` and ``ABLForcing``
   source terms act in the presence of water. However, this parameter is not needed if the ``OceanWaves`` physics
   module is being used because; the ``water_level`` value is automatically populated in that case from the ``OceanWaves``
   setup parameters.

Notes
-----

Terrain-blanked cells
^^^^^^^^^^^^^^^^^^^^^

When the ``terrain_blank`` int field is present (i.e., when a physics module such as
:doc:`TerrainDrag <inputs_terraindrag>` or ``ChannelBuilder`` is active), the volume-of-fluid
post-solve step extrapolates the vof field into the cells that are inside the terrain body
(``terrain_blank = 1``). This keeps the interface within the solid from drifting away from the
surrounding flow, which would otherwise lead to spurious density gradients and forcing near the
terrain surface.

The extrapolation is performed after the vof solve and before the density is computed from vof:
- Blanked cells with at least one lateral neighbor (in x or y) that is not blanked take the average
   vof of those unblanked neighbors. Valid cells on the x or y domain boundary may be used as
   sources, but ghost cells outside the domain are excluded.
- Blanked cells that have no non-blanked lateral neighbors take the average of the laterally
  neighboring blanked cells that have already been assigned a value. This is repeated in sweeps,
  so values propagate one layer deeper into the terrain per sweep until every reachable blanked
  cell has been assigned.
- Only lateral (x and y) neighbors are used, so each horizontal plane of cells is filled
  independently of the planes above and below it.
- Blanked cells that are not reachable from any non-blanked cell keep their existing vof value.

This step is a no-op when no ``terrain_blank`` field exists, and it does not modify the vof values
of non-blanked cells.
