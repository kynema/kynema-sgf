.. _inputs_forestdrag:

Section: ForestDrag
~~~~~~~~~~~~~~~~~~~

These parameters are active when ``ForestDrag`` is included in
:input_param:`incflo.physics`. Two forest representations are supported:

1. A legacy cylinder-based forest file, controlled by ``ForestDrag.forest_file``.
2. A point-cloud forest description, controlled by ``ForestDrag.point_cloud_files``.

The two approaches are mutually exclusive and cannot be used together in the same run.

.. input_param:: ForestDrag.forest_file

   **type:** String, optional, default = ``forest.amrwind``

   Input file for the legacy cylinder-based forest representation. If
   ``ForestDrag`` is enabled and
   :input_param:`ForestDrag.point_cloud_files` is not provided,
   Kynema-SGF reads this file.

   Each non-comment line of the file must contain eight values,

   ``TreeType xc yc height diameter cd lai laimax``

   where ``TreeType`` selects the vertical LAD model, ``xc`` and ``yc`` are
   the forest-center coordinates, ``height`` and ``diameter`` define the
   cylindrical footprint, ``cd`` is the drag coefficient, ``lai`` is the leaf
   area index, and ``laimax`` gives the normalized height of the LAD maximum
   for the type-2 profile.

   The supported legacy LAD models are the same as those described in the
   theory documentation: a uniform profile and an analytical vertically varying
   profile.

.. input_param:: ForestDrag.point_cloud_files

   **type:** List of strings, optional

   List of point-cloud files, with one file per forest patch. When this
   parameter is present, the legacy cylinder-based input
   :input_param:`ForestDrag.forest_file` must not be specified.

   Each non-comment line of a point-cloud file must contain four values,

   ``x y z lad``

   where ``x``, ``y``, and ``z`` are the point coordinates and ``lad`` is the
   leaf area density sampled at that point.

   Blank lines are ignored. Text following ``#`` on a line is treated as a
   comment.

.. input_param:: ForestDrag.coefficients_of_drag

   **type:** List of reals, mandatory when
   :input_param:`ForestDrag.point_cloud_files` is used

   Drag coefficient for each point-cloud forest. The number of entries in this
   list must match the number of files listed in
   :input_param:`ForestDrag.point_cloud_files`.

.. input_param:: ForestDrag.point_neighbors

   **type:** Integer, optional, default = 4

   Number of nearest point-cloud samples used to reconstruct the local LAD at a
   cell center. Values are internally clamped to the range 1 through 8.

.. input_param:: ForestDrag.point_interp_eps

   **type:** Real, optional, default = 1.0e-12

   Regularization parameter used in the inverse-distance weighting for the
   point-cloud forest model. This value avoids singular weights when a cell
   center is extremely close to a point-cloud sample.