.. _inputs_damping_layer:

Section: DampingLayer
~~~~~~~~~~~~~~~~~~~~~

The DampingLayer physics module creates one or more spatially varying damping
coefficient fields near selected domain boundaries. The DampingLayerSource
source term then uses those coefficients to relax solved fields toward user-
defined targets.

Activate DampingLayer by including it in :input_param:`incflo.physics`.

Activate DampingLayerSource for each equation you want to damp, for example:

- ``ICNS.source_terms = DampingLayerSource``
- ``temperature.source_terms = DampingLayerSource``
- ``TKE.source_terms = DampingLayerSource``

Additional scalar equations that expose source-term controls can also use
``DampingLayerSource``.

.. input_param:: DampingLayer.fields

   **type:** List of strings, mandatory when ``DampingLayer`` is active

   Field labels that will have damping-layer coefficient fields created.
   Typical values include ``velocity``, ``temperature``, ``density``, ``tke``,
   and ``sdr``.

For each entry in :input_param:`DampingLayer.fields`, parameters are configured
per boundary using:

``DampingLayer.<field>.<boundary>.*`` where ``boundary`` is one of
``xlo``, ``xhi``, ``ylo``, ``yhi``, ``zlo``, or ``zhi``.

.. input_param:: DampingLayer.<field>.<boundary>.thickness

   **type:** Real, optional

   Thickness of the damping region for this boundary. A value ``> 0`` enables
   damping on that boundary for this field.

.. input_param:: DampingLayer.<field>.<boundary>.blending_fraction

   **type:** Real, optional, default = 0.0

   Fraction of ``thickness`` used as a ramp from full damping to zero damping.
   The remaining fraction applies full damping.

.. input_param:: DampingLayer.<field>.<boundary>.blending_function_type

   **type:** String, optional, default = cosine

   Shape of the ramp region. Supported values are ``linear``, ``quadratic``,
   ``exponential``, and ``cosine``.

.. input_param:: DampingLayer.<field>.<boundary>.minimum_height

   **type:** Real, optional (only for ``x*`` and ``y*`` boundaries)

   Adds a vertical limiter so horizontal boundary damping is only applied above
   this height. This parameter is not valid on ``zlo`` or ``zhi``.

.. input_param:: DampingLayer.<field>.<boundary>.vertical_blending_thickness

   **type:** Real, mandatory when ``minimum_height`` is provided

   Vertical blending thickness used by the vertical limiter activated through
   ``minimum_height``. This parameter is only valid on ``x*`` and ``y*``
   boundaries.

.. input_param:: DampingLayer.<field>.<boundary>.vertical_blending_function_type

   **type:** String, optional, default = cosine

   Shape for blending in the vertical direction. Supported values are ``linear``,
   ``quadratic``, ``exponential``, and ``cosine``. This parameter is only
   valid on ``x*`` and ``y*`` boundaries.

DampingLayerSource target parameters are configured with the same namespace:

``DampingLayer.<field>.<boundary>.*``

.. input_param:: DampingLayer.<field>.<boundary>.target_type

   **type:** String, mandatory when damping is active on a boundary

   Target mode for the relaxed value. Supported values are ``constant``,
   ``profile``, ``function``, and ``field``.

.. input_param:: DampingLayer.<field>.<boundary>.target_value

   **type:** List of Real, mandatory for ``target_type = constant``

   Constant target value per component. The number of values must match the
   number of components in the damped field.

.. input_param:: DampingLayer.<field>.<boundary>.target_profile_heights

   **type:** List of Real, mandatory for ``target_type = profile``

   Monotone height coordinates used for profile interpolation.

.. input_param:: DampingLayer.<field>.<boundary>.target_profile_values

   **type:** List of Real, optional for ``target_type = profile``

   Convenience form for scalar or x-component profile values.

.. input_param:: DampingLayer.<field>.<boundary>.target_profile_values_x

   **type:** List of Real, optional for ``target_type = profile``

   x-component profile values.

.. input_param:: DampingLayer.<field>.<boundary>.target_profile_values_y

   **type:** List of Real, optional for ``target_type = profile``

   y-component profile values.

.. input_param:: DampingLayer.<field>.<boundary>.target_profile_values_z

   **type:** List of Real, optional for ``target_type = profile``

   z-component profile values.

.. input_param:: DampingLayer.<field>.<boundary>.target_function

   **type:** String expression, mandatory for ``target_type = function``

   Expression parsed at runtime. Available variables are ``t``, ``x``, ``y``,
   ``z``, and ``n`` (component index).

.. input_param:: DampingLayer.<field>.<boundary>.target_field

   **type:** String, mandatory for ``target_type = field``

   Name of another field used as the target value. This could work with a field
   created by a utility or post-processing routine, such as an averaging field,
   or a field created by a dedicated physics module.

.. input_param:: DampingLayer.<field>.<boundary>.damped_components

   **type:** List of Integer, optional

   Per-component on/off mask for damping. The list length must match the
   number of components in the damped field, and each entry must be either
   ``0`` (off) or ``1`` (on). If omitted, all components default to ``1``.

Example
^^^^^^^

.. code-block:: console

   incflo.physics = ABL DampingLayer
   DampingLayer.fields = velocity temperature

   ICNS.source_terms = BoussinesqBuoyancy DampingLayerSource
   temperature.source_terms = DampingLayerSource

   DampingLayer.velocity.xlo.thickness = 150.0
   DampingLayer.velocity.xlo.blending_fraction = 0.25
   DampingLayer.velocity.xlo.blending_function_type = linear
   DampingLayer.velocity.xlo.minimum_height = 100.0
   DampingLayer.velocity.xlo.vertical_blending_thickness = 80.0
   DampingLayer.velocity.xlo.vertical_blending_function_type = cosine
   DampingLayer.velocity.xlo.target_type = constant
   DampingLayer.velocity.xlo.target_value = 8.0 0.0 0.0
   DampingLayer.velocity.xlo.damped_components = 1 0 0

   DampingLayer.velocity.ylo.thickness = 120.0
   DampingLayer.velocity.ylo.blending_fraction = 0.5
   DampingLayer.velocity.ylo.blending_function_type = quadratic
   DampingLayer.velocity.ylo.target_type = profile
   DampingLayer.velocity.ylo.target_profile_heights = 0.0 500.0 1000.0
   DampingLayer.velocity.ylo.target_profile_values_x = 6.0 7.0 8.0
   DampingLayer.velocity.ylo.target_profile_values_y = 0.0 0.5 1.0
   DampingLayer.velocity.ylo.target_profile_values_z = 0.0 0.0 0.0

   DampingLayer.temperature.zhi.thickness = 200.0
   DampingLayer.temperature.zhi.blending_fraction = 0.4
   DampingLayer.temperature.zhi.blending_function_type = cosine
   DampingLayer.temperature.zhi.target_type = function
   DampingLayer.temperature.zhi.target_function = "300.0 + 0.01*z"
