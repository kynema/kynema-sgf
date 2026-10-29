.. _inputs_fieldnorms:

Section: FieldNorms
~~~~~~~~~~~~~~~~~~~

This section controls field norms post-processing. Field norms
output the global norm of every plot variable to a text file, writing
a new line for each output time.
The prefix is the label set in ``incflo.post_processing``. For example
``incflo.post_processing = fieldnorms``. Inputs listed here only address
options specific to field norms; the inputs controlling the output timing
(which are shared by other post-processing types) are listed in the
[:ref:`post-processing section <inputs_post_processing>`].

.. input_param:: fieldnorms.type

   **type:** String, mandatory

   To use field norms output specify with keyword ``FieldNorms``
   
.. input_param:: fieldnorms.mask_redundant_grids

   **type:** Boolean, optional, default = true

   Setting this option to true ensures that the field norm calculation
   does not include cells that are covered by finer cells from a different
   mesh level. This also ensures that no volume is double counted. Setting
   it to false makes sure that every cell, covered or not, is included.

.. input_param:: fieldnorms.use_vector_magnitude

   **type:** Boolean, optional, default = false

   By default, each component of vector output fields is considered separately.
   Setting this option to true ensures that the field norm calculation
   considers the local vector magnitude instead of separating each vector component.

.. input_param:: fieldnorms.norm_type

   **type:** String, optional, default = ``2``

   By default, this post-processing routine calculates the volume-weighted
   L2 norm, corresponding to norm type 2. Other options include the L1 norm,
   designated by ``1``, and the L-infinity norm, designated by ``infinity``.