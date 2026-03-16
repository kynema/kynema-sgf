.. _inputs_subvolume:
  
Section: Subvolume
~~~~~~~~~~~~~~~~~~

This section controls subvolume post-processing. Subvolume outputs chunks
of data directly from the computational mesh.
The prefix is the label set in ``incflo.post_processing``. For example
``incflo.post_processing = subvol1``


.. input_param:: subvol1.type

   **type:** String, mandatory

   To use subvolume post-processing, specify with keyword ``Subvolume``

.. input_param:: subvol1.labels

   **type:** List of strings, mandatory

   Similar to the Sampling utility, multiple subvolumes can be defined within a
   single Subvolume instance, where the top-level label is used to define the fields
   and output parameters and the bottom-level labels are used to define the type
   of subvolume and its spatial parameters. Below, the label ``chunk1`` will be used
   where applicable.

.. input_param:: subvol1.fields

   **type:** List of strings, mandatory

   Specify which field arrays should be output within the selected subvolume 

.. input_param:: subvol1.int_fields

   **type:** List of strings, optional, default is empty

   Specify which integer field arrays should be output within the selected subvolume

.. input_param:: subvol1.derived_fields

   **type:** List of strings, optional, default is empty

   Specify which derived field arrays should be output within the selected subvolume

.. input_param:: subvol1.output_rename

   **type:** String, optional

   If desired, use a different name for the top-level label when writing. In this example,
   this would replace the ``subvol1`` label in the naming of the output directories. This
   option is primarily intended to enable different subvolumes to have the same top-level
   names despite needing different top-level parameters. In that case, it is the responsibility
   of the user to ensure the bottom-level names are still different, as this is not checked by the code.

.. input_param:: subvol1.chunk1.type

   **type:** String, optional, default = Rectangular

   This specifies the type of subvolume to be used. At the moment, Rectangular is the only available type.
   Rectangular subvolumes are defined using an origin, number of cells, and cell size.

.. input_param:: subvol1.chunk1.origin

   **type:** Vector<Real>, mandatory

   Starting point in three-dimensional space to define the rectangular subvolume. This needs
   to be the lower left corner of a mesh cell for the subvolume data extraction to work properly.

.. input_param:: subvol1.chunk1.num_points

   **type:** Vector<Int>, mandatory

   Number of points (more accurately, cells) in each direction to include in the rectangular subvolume.

.. input_param:: subvol1.chunk1.dx_vec

   **type:** Vector<Real>, mandatory

   Cell size, in each direction, for the rectangular subvolume. This is used to determine which mesh level
   should be used for the subvolume, and, as a result, this input argument must correspond to the resolution
   of one of the mesh levels in the spatial extent of the subvolume. If a single cell size is specified through
   the :input_param:`subvol1.chunk1.dx` argument below, then this argument (:input_param:`subvol1.chunk1.dx_vec`) is not required.

.. input_param:: subvol1.chunk1.dx

   **type:** Real, optional

   Cell size for the rectangular subvolume. This optional argument can be used in place of
   :input_param:`subvol1.chunk1.dx_vec` to specify a single cell size, representing the
   cell size in all three directions.

.. input_param:: subvol1.chunk1.chunk_size_vec

   **type:** Vector<Int>, optional

   Optional argument for more control over how the data is partitioned within a subvolume. By default,
   the chunk size will be the max grid size in each spatial direction.
