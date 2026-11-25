.. _inputs_geometry:

Section: geometry
~~~~~~~~~~~~~~~~~

This section deals with inputs related to the problem domain.

.. input_param:: geometry.prob_lo

   **type:** List of 3 real numbers, mandatory

   The coordinates of *lower corner* of the computational domain bounding box.

.. input_param:: geometry.prob_hi

   **type:** List of 3 real numbers, mandatory

   The coordinates of the *upper corner* of the computational domain bounding box.
   
.. input_param:: geometry.prob_hi_physical

  **type:** List of 3 real numbers, optional, default = ``geometry.prob_hi``

  The coordinates of the *upper corner* of the computational domain bounding box for a
  Stretched Mesh.

.. input_param:: geometry.is_periodic

   **type:** List of 3 integers, mandatory

   Flags indicating whether the flow is periodic in the ``x``, ``y``, or ``z``
   directions respectively.
   
.. input_param:: geometry.mesh_mapping

 **type:** String, optional, default = ``ConstantMap``

 Define the style of mapping between the stretched coordinates and the uniform coordinates. The default map is a constant scaling map with 
 ``ConstantMap.scaling_factor = 1.0 1.0 1.0``.
