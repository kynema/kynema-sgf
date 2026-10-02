.. _inputs_terraindrag:

Section: TerrainDrag
~~~~~~~~~~~~~~~~~~~~

These parameters are active when ``TerrainDrag`` is included in
:input_param:`incflo.physics`.

TerrainDrag has two operation modes:

1. Terrain-file mode (default), where terrain and roughness are read from
   files.
2. Ocean-wave mode, activated automatically when ``OceanWaves`` physics is
   active and the ``vof`` field is not present. In this mode, terrain blanking
   and drag markers are constructed from wave fields.

.. input_param:: TerrainDrag.terrain_file

   **type:** String, optional, default = ``terrain.amrwind``

   Input file for terrain height data in terrain-file mode.

.. input_param:: TerrainDrag.roughness_file

   **type:** String, optional, default = ``terrain.roughness``

   Input file for roughness-length data in terrain-file mode.

   If this file is missing or cannot be opened, the roughness is set from
   :input_param:`TerrainDrag.uniform_roughness`.

.. input_param:: TerrainDrag.uniform_roughness

   **type:** Real, optional, default = 0.1

   Uniform roughness length used when
   :input_param:`TerrainDrag.roughness_file` is not available.

   If both ``uniform_roughness`` and ``roughness_file`` are provided, values
   from the roughness file are used where available.

.. input_param:: TerrainDrag.damp_east_slope

   **type:** Real, optional, default = 0.0

   Width of the east (+x) damping ramp region.

.. input_param:: TerrainDrag.damp_east_full

   **type:** Real, optional, default = 0.0

   Width of the east (+x) full-strength damping region adjacent to the boundary.

.. input_param:: TerrainDrag.damp_west_slope

   **type:** Real, optional, default = 0.0

   Width of the west (-x) damping ramp region.

.. input_param:: TerrainDrag.damp_west_full

   **type:** Real, optional, default = 0.0

   Width of the west (-x) full-strength damping region adjacent to the boundary.

.. input_param:: TerrainDrag.damp_north_slope

   **type:** Real, optional, default = 0.0

   Width of the north (+y) damping ramp region.

.. input_param:: TerrainDrag.damp_north_full

   **type:** Real, optional, default = 0.0

   Width of the north (+y) full-strength damping region adjacent to the boundary.

.. input_param:: TerrainDrag.damp_south_slope

   **type:** Real, optional, default = 0.0

   Width of the south (-y) damping ramp region.

.. input_param:: TerrainDrag.damp_south_full

   **type:** Real, optional, default = 0.0

   Width of the south (-y) full-strength damping region adjacent to the boundary.

.. input_param:: TerrainDrag.horizontal_time_scale

   **type:** Real, optional, default = 20.0

   Time scale used to convert the assembled damping coefficient into
   ``terrain_damping`` forcing (larger values reduce damping strength).

.. input_param:: TerrainDrag.horizontal_abl_height

   **type:** Real, optional, default = 1.0e15

   Height at which horizontal damping begins to turn on.

.. input_param:: TerrainDrag.horizontal_slope_end

   **type:** Real, optional, default = 1.0e30

   Height where the horizontal damping vertical multiplier reaches full value.

.. input_param:: TerrainDrag.vertical_slope

   **type:** Real, optional, default = 0.0

   Height where additional full-domain vertical damping begins to ramp on.

.. input_param:: TerrainDrag.vertical_full

   **type:** Real, optional, default = 1.0e30

   Height where additional full-domain vertical damping reaches full strength.

Notes
-----

- The terrain drag marker field (``terrain_drag``) is derived from transitions
  in ``terrain_blank`` and is used by :doc:`inputs_Momentum_Sources` when
  ``DragForcing`` is active.
- In ocean-wave mode, TerrainDrag derives ``terrain_blank`` and
  ``terrain_height`` from wave fields and sets a uniform low roughness
  (``terrainz0 = 1.0e-4``).
- In ocean-wave mode, terrain-file and roughness-file parameters are not used.
