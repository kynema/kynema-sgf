.. _inputs_transport:

Section: transport
~~~~~~~~~~~~~~~~~~

This section is for setting thermal and momentum diffusivity coefficients.

.. input_param:: transport.model

   **type:** String, optional, default = ConstTransport

   The implemented transport model are ConstTransport and TwoPhaseTransport.
   
.. input_param:: transport.viscosity

   **type:** Real, optional, default = 1.0e-5

   Sets the dynamic viscosity
   
.. input_param:: transport.laminar_prandtl 

   **type:** Real, optional, default = 1.0

   Sets the laminar Prandtl number.
   
.. input_param:: transport.turbulent_prandtl 

   **type:** Real, optional, default = 1.0

   Sets the turbulent Prandtl number.
   
.. input_param:: transport.reference_temperature

   **type:** Real, mandatory
   
   Reference temperature :math:`\theta_\mathrm{ref}` in Kelvin.
   Values of the temperature field that are less than or greater than this value will 
   cause a buoyancy force in the direction of the gravity vector.
   
.. input_param:: transport.thermal_expansion_coefficient

   **type:** Real, optional, default :math:`\beta = 1 / \theta_\mathrm{ref}`
   
   Thermal expansion coefficient, if not specified this value is set to the inverse of the
   :input_param:`transport.reference_temperature` value.
   
