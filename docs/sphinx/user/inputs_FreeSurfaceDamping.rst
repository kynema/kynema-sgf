.. _inputs_FreeSurfaceDamping:

Section: FreeSurfaceDamping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section controls the free-surface damping physics module. The module is
intended for multiphase simulations and applies relaxation near the domain
boundaries to reduce wave reflections and smooth the free surface.

Example:

.. code-block:: console

   incflo.physics = MultiPhase FreeSurfaceDamping
   FreeSurfaceDamping.vertical_velocity_damping = true
   FreeSurfaceDamping.volume_fraction_damping = true
   FreeSurfaceDamping.global_damping = false
   FreeSurfaceDamping.length_xlo = 5.0
   FreeSurfaceDamping.length_xhi = 10.0
   FreeSurfaceDamping.length_ylo = 5.0
   FreeSurfaceDamping.length_yhi = 5.0

.. input_param:: FreeSurfaceDamping.vertical_velocity_damping

   **type:** Boolean, optional, default = true

   When enabled, the module damps the vertical velocity component toward zero
   in the active damping region.

.. input_param:: FreeSurfaceDamping.volume_fraction_damping

   **type:** Boolean, optional, default = true

   When enabled, the module smooths the volume fraction field using a local
   horizontal 3x3 average before relaxing it toward the current state. If the
   volume fraction changes, the density field is updated consistently.

.. input_param:: FreeSurfaceDamping.global_damping

   **type:** Boolean, optional, default = true

   When enabled, damping is applied throughout the domain. When disabled,
   damping is restricted to boundary zones whose widths are specified by the
   ``length_x*`` and ``length_y*`` parameters. Within those boundary zones,
   the damping strength is ramped up from zero at the edge of the damping zone
   to full strength at the domain boundary.

.. input_param:: FreeSurfaceDamping.time_scale_fraction

   **type:** Real, optional, default = 0.5

   Fraction of the computed damping update to apply on each call to the
   method. Smaller values make the damping act more gradually. This parameter
   should not exceed 1.

.. input_param:: FreeSurfaceDamping.length_xlo

   **type:** Real, mandatory when ``global_damping = false``

   Physical length of the damping zone near the low-x boundary.

.. input_param:: FreeSurfaceDamping.length_xhi

   **type:** Real, mandatory when ``global_damping = false``

   Physical length of the damping zone near the high-x boundary.

.. input_param:: FreeSurfaceDamping.length_ylo

   **type:** Real, mandatory when ``global_damping = false``

   Physical length of the damping zone near the low-y boundary.

.. input_param:: FreeSurfaceDamping.length_yhi

   **type:** Real, mandatory when ``global_damping = false``

   Physical length of the damping zone near the high-y boundary.
