.. _inputs_tagging:

Section: tagging
~~~~~~~~~~~~~~~~

This section manages the various mesh refinement criteria that can be used to
activate either static or adaptive mesh refinement during simulations. The
parameters are read from the prefix ``tagging`` and can contain different types
of tagging logic. Note that this section is only active if
:input_param:`amr.max_level` is greater than zero. Regridding interval is controlled by :input_param:`time.regrid_interval` .

Example::

  tagging.labels = s1 f1 g1
  tagging.s1.type = CartBoxRefinement
  tagging.s1.static_refinement_def = static_box.txt

  tagging.f1.type = FieldRefinement
  tagging.f1.field_name = density
  tagging.f1.grad_error = 0.1. 0.1 0.1

  tagging.g1.type = GeometryRefinement
  tagging.g1.shapes = c1 b1

  tagging.g1.c1.type = cylinder
  tagging.g1.c1.start = 500.0 500.0 250.0
  tagging.g1.c1.end = 500.0 500.0 750.0
  tagging.g1.c1.outer_radius = 300.0
  tagging.g1.c1.inner_radius = 275.0

  tagging.g1.b1.type = box
  tagging.g1.b1.origin = 300.0 150.0 250.0
  tagging.g1.b1.xaxis =  450.0 600.0 0.0
  tagging.g1.b1.yaxis =  -150.0 100.0 0.0
  tagging.g1.b1.zaxis = 0.0 0.0 500.0

Each section must contain the keyword ``type`` that is one of the refinement types:

========================== ===================================================================
``CartBoxRefinement``      Nested refinement using Cartesian boxes
``FieldRefinement``        Refinement based on error metric for field or its gradient
``OversetRefinement``      Refinement around fringe/field interface
``GeometryRefinement``     Refinement using geometric shapes
``QCriterionRefinement``   Refinement using Q-Criterion
``VorticityMagRefinement`` Refinement using vorticity
========================== ===================================================================

.. input_param:: tagging.labels

   **type:** List of one or more names

   Labels indicate a list of prefixes for different types of refinement criteria
   active during the simulation.

The parameters for the subsections are determined by the type of refinement being performed.

.. _inputs_static_refinement:

Refinement using Cartesian boxes
````````````````````````````````

``CartBoxRefinement`` allows refining boxes (aligned with the principal axes).

Example::

   tagging.labels = static
   tagging.static.type = CartBoxRefinement
   tagging.static.static_refinement_def = static_box.txt

.. input_param:: tagging.CartBoxRefinement.static_refinement_def

   **type:** String, required

   The text file that contains a list of bounding boxes used to perform
   refinement at various levels.

Refinement using geometry
`````````````````````````

This section controls refinement using pre-defined geometric shapes. Currently,
two options are supported: 1. ``box`` -- refines the region inside a hexahedral
block, 2. ``cylinder`` -- refines the region inside a cylindrical block, and
3. ``udf`` -- refines along an analytical function depending on time and spatial coordinate.
With the exception of the ``udf`` able to prescribe temporal variations, geometry-based refinement is static.

.. input_param:: tagging.GeometryRefinement.shapes

   **type:** List of strings, required

   Names of the input subsections that define specific geometries for refinement.

.. input_param:: tagging.GeometryRefinement.level

   **type:**  Integer, optional, default: -1

   If ``level`` is provided and is greater than or equal to 0, then the
   refinement based on geometries defined for this section is only performed at
   that level.

.. input_param:: tagging.GeometryRefinement.min_level

   **type:**  Integer, optional, default: 0

   If ``level`` is not specified, then this option specifies the minimum level
   where this refinement is active.

.. input_param:: tagging.GeometryRefinement.max_level

   **type:**  Integer, optional, default: ``mesh.maxLevel()``

   If ``level`` is not specified, then this option specifies the maximum level
   where this refinement is active.

Note that the specification of ``level`` overrides, ``min_level`` and
``max_level`` specifications. This can be used to control the different levels
where refinement regions are active.

Example::

  tagging.g1.type = GeometryRefinement
  tagging.g1.shapes = b1 b2
  tagging.g1.level = 0
  tagging.g1.b1.type = box
  tagging.g1.b1.origin = 300.0 150.0 250.0
  tagging.g1.b1.xaxis =  450.0 600.0 0.0
  tagging.g1.b1.yaxis =  -150.0 100.0 0.0
  tagging.g1.b1.zaxis = 0.0 0.0 500.0
  tagging.g1.b2.type = box
  tagging.g1.b2.origin = 600.0 350.0 250.0
  tagging.g1.b2.xaxis =  50.0 30.0 0.0
  tagging.g1.b2.yaxis =  -50.0 60.0 0.0
  tagging.g1.b2.zaxis = 0.0 0.0 500.0

  tagging.g2.type = GeometryRefinement
  tagging.g2.shapes = c1
  tagging.g2.level = 1
  tagging.g2.c1.type = cylinder
  tagging.g2.c1.start = 500.0 500.0 250.0
  tagging.g2.c1.end = 500.0 500.0 750.0
  tagging.g2.c1.outer_radius = 300.0
  tagging.g2.c1.inner_radius = 275.0

  tagging.g3.type = GeometryRefinement
  tagging.g3.shapes = udf0
  tagging.g3.level = 0
  tagging.g3.udf0.type = udf
  tagging.g3.udf0.udf = "if(x > 500, 1.0, 0.0)"
  tagging.g3.udf0.box_lo = 0.0 0.0 0.0
  tagging.g3.udf0.box_hi = 1000.0 1000.0 500.0


This example defines three different refinement definitions acting on
levels 0, 1 and 0 respectively. The first refinement at level 0
(``g1``) contains two box regions, whereas the refinement at level 1
(``g2``) only contains one cylinder definition. The second refinement
at level 0 (``g3``) uses a user defined function (UDF) to specify a
geometrical refinement region.

**Refinement using hexahedral block definitions**

To perform ``box`` refinement, the user specifies the ``origin`` of the box and
three vectors: ``xaxis, yaxis, zaxis`` that defines the directions and the
extents of the hexahedral block. Denoting :math:`\mathbf{O}` as origin vector
and :math:`\mathbf{x}`, :math:`\mathbf{y}` and :math:`\mathbf{z}` as the three
vectors given by the user, the position vectors of the eight corners of the
hexahedral box are given by

.. math::

   \mathbf{x}_0 &= \mathbf{O} && \mathbf{x}_4 &= \mathbf{O} + \mathbf{z} \\
   \mathbf{x}_1 &= \mathbf{O} + \mathbf{x} && \mathbf{x}_5 &= \mathbf{O} + \mathbf{z} + \mathbf{x} \\
   \mathbf{x}_2 &= \mathbf{O} + \mathbf{x} + \mathbf{y} \qquad && \mathbf{x}_6 &= \mathbf{O} + \mathbf{z} + \mathbf{x} + \mathbf{y} \\
   \mathbf{x}_3 &= \mathbf{O} + \mathbf{y} && \mathbf{x}_7 &= \mathbf{O} + \mathbf{z} + \mathbf{y} \\



**Refinement using cylindrical block definitions**

The axis and the extents along the axis are defined by two position vectors
``start`` and ``end``. The radial extent is specified by ``outer_radius``. An
optional ``inner_radius`` can be specified to restrict tagging to an annulus
between the inner and outer radii.

**Refinement using a user specified function**

The UDF uses AMReX's Parser class to parse a string in the input file,
``udf`` key in the input file. For more information on the AMReX
Parser and supported functions, consult the `AMReX Parser
documentation
<https://amrex-codes.github.io/amrex/docs_html/Basics.html#parser>`_.
The UDF can be (almost) arbitrarily complex, leveraging standard
functions, boolean operators, local variables, etc. The following
assumptions are imposed on the UDF: 1. it must evaluate to 0 (tagging
off) or 1 (tagging on) (the return value of the UDF is cast to a
boolean and that is used to set the tags), 2. the allowed variables
are ``t``, ``x``, ``y``, and ``z``.

.. input_param:: tagging.GeometryRefinement.UDFRefiner.udf

   **type:** String, required

   String specifying the user defined function.

.. input_param:: tagging.GeometryRefinement.UDFRefiner.box_lo

   **type:** Vector<Real>, optional

   List of the low corner values for a bounding box where the tagging
   will be active. By default the bounding box will span the entire domain.

.. input_param:: tagging.GeometryRefinement.UDFRefiner.box_hi

   **type:** Vector<Real>, optional

   List of the high corner values for a bounding box where the tagging
   will be active. By default the bounding box will span the entire domain.

Actuator-following refinement
`````````````````````````````

``ActuatorRefinement`` maintains rotor-aligned cylindrical refinement regions
around selected actuator models. A selected ``ActuatorSector`` contributes one
cylinder. A selected ``Drone`` contributes one cylinder for each child rotor;
overlapping cylinders are combined by the normal AMR tagging operation.

The cylinder center and axis follow the prescribed rotor-hub position and
rotor normal. ``forward`` always means the positive rotor-normal direction and
``backward`` means the negative direction. These names are geometric and do
not imply upstream or downstream.

For example, the following applies one configuration to every rotor belonging
to drones ``Q1`` and ``Q2``::

  tagging.labels = rotor_tracking
  tagging.rotor_tracking.type = ActuatorRefinement
  tagging.rotor_tracking.actuator_labels = Q1 Q2
  tagging.rotor_tracking.min_level = 0
  tagging.rotor_tracking.max_level = 1
  tagging.rotor_tracking.radial_padding_epsilon = 3.0
  tagging.rotor_tracking.forward_padding_diameter = 1.0
  tagging.rotor_tracking.backward_padding_diameter = 0.5

``min_level`` and ``max_level`` identify levels on which cells are tagged, so
``max_level`` must be less than ``amr.max_level``.

.. input_param:: tagging.ActuatorRefinement.actuator_labels

   **type:** List of strings, mandatory

   Top-level labels from ``Actuator.labels``. Each label must identify an
   ``ActuatorSector`` or ``Drone``. Selecting a drone automatically selects all
   child rotors.

.. input_param:: tagging.ActuatorRefinement.min_level

   **type:** Integer, optional, default = 0

   Lowest mesh level on which this criterion tags cells.

.. input_param:: tagging.ActuatorRefinement.max_level

   **type:** Integer, optional, default = ``amr.max_level - 1``

   Highest mesh level on which this criterion tags cells.

.. input_param:: tagging.ActuatorRefinement.radial_padding_epsilon

   **type:** Real number, optional, default = 3.0

   Additional cylinder radius as a multiple of the rotor's maximum Gaussian
   width. The total radius also includes the physical rotor radius.

.. input_param:: tagging.ActuatorRefinement.radial_padding

   **type:** Real number, optional, default = 0.0

   Additional cylinder radius in meters.

.. input_param:: tagging.ActuatorRefinement.axial_padding_epsilon

   **type:** Real number, optional, default = 3.0

   Cylinder extent on both sides of the rotor plane as a multiple of maximum
   Gaussian width. A side-specific epsilon input overrides this value on that
   side.

.. input_param:: tagging.ActuatorRefinement.axial_padding_diameter

   **type:** Real number, optional

   Cylinder extent on both sides of the rotor plane in rotor-diameter units. A
   side-specific diameter input overrides this value on that side.

.. input_param:: tagging.ActuatorRefinement.axial_padding

   **type:** Real number, optional, default = 0.0

   Cylinder extent in meters on both sides of the rotor plane. A side-specific
   absolute input overrides this value on that side.

.. input_param:: tagging.ActuatorRefinement.forward_padding_epsilon

   **type:** Real number, optional

   Forward extent as a multiple of maximum Gaussian width. If omitted,
   ``axial_padding_epsilon`` is used.

.. input_param:: tagging.ActuatorRefinement.backward_padding_epsilon

   **type:** Real number, optional

   Backward extent as a multiple of maximum Gaussian width. If omitted,
   ``axial_padding_epsilon`` is used.

.. input_param:: tagging.ActuatorRefinement.forward_padding_diameter

   **type:** Real number, optional

   Forward extent in rotor-diameter units. If omitted,
   ``axial_padding_diameter`` is used.

.. input_param:: tagging.ActuatorRefinement.backward_padding_diameter

   **type:** Real number, optional

   Backward extent in rotor-diameter units. If omitted,
   ``axial_padding_diameter`` is used.

.. input_param:: tagging.ActuatorRefinement.forward_padding

   **type:** Real number, optional

   Forward extent in meters. If omitted, ``axial_padding`` is used.

.. input_param:: tagging.ActuatorRefinement.backward_padding

   **type:** Real number, optional

   Backward extent in meters. If omitted, ``axial_padding`` is used.

Axial epsilon, diameter, and absolute inputs are mutually exclusive. Specify
only one unit family for a refinement criterion. Within that family, a shared
``axial_padding_*`` value can be combined with forward and backward overrides.
If no axial input is supplied, both sides default to ``3 epsilon``.


Refinement using field error criteria
`````````````````````````````````````

Example::

  tagging.f1.type = FieldRefinement
  tagging.f1.field_name = density
  tagging.f1.grad_error = 0.1 0.1 0.1
  tagging.f1.box_lo = 10.0 10.0 10.0
  tagging.f1.box_hi = 20.0 20.0 20.0

.. input_param:: tagging.FieldRefinement.field_name

   **type:** String, required

   The name of the field used to tag cells

.. input_param:: tagging.FieldRefinement.field_error

   **type:** Vector<Real>, optional

   List of field error values at each level. The user must specify a value for
   each level desired.

.. input_param:: tagging.FieldRefinement.grad_error

   **type:** Vector<Real>, optional

   List of gradient error values at each level. The user must specify a value for
   each level desired.

.. input_param:: tagging.FieldRefinement.box_lo

   **type:** Vector<Real>, optional

   List of the low corner values for a bounding box where the tagging
   will be active. By default the bounding box will span the entire domain.

.. input_param:: tagging.FieldRefinement.box_hi

   **type:** Vector<Real>, optional

   List of the high corner values for a bounding box where the tagging
   will be active. By default the bounding box will span the entire domain.


Refinement using Q-Criterion
`````````````````````````````````````

Example::

  tagging.qc1.type = QCriterionRefinement
  tagging.qc1.nondim = false
  tagging.qc1.values = 10.0 20.0 20.0

.. input_param:: tagging.QCriterionRefinement.nondim

   **type:** Boolean, optional, default = true

   Boolean determining if the dimensional or non-dimensional form 
   of Q-criterion should be used. Dimensional version may require 
   modifying values depending on physical scales. For the non-dimensional 
   form positive thresholds indicate regions where the rotational strength is 
   larger than the shear rate strength. A threshold of unity indicates 
   that the rotational strength is equal to the background shear strength. 
   
.. input_param:: tagging.QCriterionRefinement.values

   **type:** Vector<Real>, optional

   List of Q-criterion values at each level.
   If the absolute value of Q-criterion exceeds this value
   the cell is tagged for refinement.
   The user must specify a value for each level desired.

Combining refinement criteria with logical operators
````````````````````````````````````````````````````

.. input_param:: tagging.<label>.operator

   **type:** String, optional, default: ``or``

   Logical operator used to combine the tags generated by this criterion with
   the existing tag state from previously applied criteria. The operators
   are applied in the order of the labels specified in :input_param:`tagging.labels`.

   Allowed values are:

================== ====================================================================
``and``            Keep cells tagged only where both previous and current are tagged
``or``             Tag cells where either previous or current is tagged
``and_not``        Keep cells tagged where previous is tagged and current is not tagged
``or_not``         Tag cells where previous is tagged or current is not tagged
================== ====================================================================

   Note: Some refinement criteria may restrict supported operators (e.g., ``CartBoxRefinement`` currently only supports ``or``).

   In this context, previous refers to the tag state resulting from all criteria applied
   before the current one.

Example::

  tagging.labels = base hole geometry

  tagging.base.type = CartBoxRefinement
  tagging.base.operator = or
  tagging.base.static_refinement_def = static_box.txt

  tagging.hole.type = GeometryRefinement
  tagging.hole.operator = and_not
  tagging.hole.shapes = c1
  tagging.hole.c1.type = cylinder
  tagging.hole.c1.start = 500.0 500.0 250.0
  tagging.hole.c1.end = 500.0 500.0 750.0
  tagging.hole.c1.outer_radius = 300.0

  tagging.geometry.type = FieldRefinement
  tagging.geometry.operator = and
  tagging.geometry.field_name = density
  tagging.geometry.grad_error = 0.1 0.1 0.1

In this example, ``base`` creates an initial tagged region.
It is always best to start with an ``or`` operator for the initial tag.
Then, ``hole`` removes cells inside the cylinder from that region (using ``and_not``).
Finally, ``geometry`` adds field-based tags into the region (using ``and``)
that has been defined by ``base`` and modified by ``hole``.
