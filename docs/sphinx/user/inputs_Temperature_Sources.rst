.. _inputs_temperature_sources:
   
Section: Temperature Sources
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
.. input_param:: temperature.source_terms

   **type:** String(s), optional
   
   Activates source terms for the energy equations. These strings can be 
   entered in any order with a space between
   each. Please consult the :doc:`../doxygen/html/index` for a
   comprehensive list of all energy source terms available. Note that the
   following input arguments specific to each source term will only be active
   if the corresponding source term (the root name) is listed in 
   :input_param:`temperature.source_terms`.

.. input_param:: DragTempForcing.drag_coefficient

   **type:** Real, optional

   This value specifies the coefficient for the forcing term in the immersed boundary forcing method. It is currently
   recommended to use the default value to avoid initial numerical stability. 

.. input_param:: DragTempForcing.bc_forcing_time_factor

   **type:** Real, optional, default = 5.0

   This value modifies the time scale of the BC forcing component of DragTempForcing relative to
   the time step size.


The following list of inputs are used with the `Temperature.source_terms = PerturbationForcing` option to add perturbation to the 
temperature field to generate flow structures for LES when the inflow data is coarse or uniform flow condition. Not 
recommended for use with RANS models. 

.. input_param:: PerturbationForcing.start

   **type:** Real, mandatory

   Start location of the perturbation box 

.. input_param:: PerturbationForcing.end

   **type:** Real, mandatory

   End location of the perturbation box 

..  input_param:: PerturbationForcing.pert_amplitude

   **type:** Real, optional 

   Amplitude of temperature perturbation 

..  input_param:: PerturbationForcing.time_steps 

   **type:** Real, optional 

   Separation time between applying perturbations. A high value may dampen the flow structures 
   and a small value may cause numerical instability. 

.. input_param:: EBDragTempForcing.drag_coefficient

   **type:** Real, optional, default = 1.0

   This value specifies the drag/heat exchange coefficient for the EBDrag temperature forcing term.

.. input_param:: EBDragTempForcing.soil_temperature

   **type:** Real, optional, default = 300.0

   This value specifies the constant soil temperature used when neither IDW nor heatflux models are enabled.

.. input_param:: EBDragTempForcing.bc_forcing_time_factor

   **type:** Real, optional, default = 5.0

   This value modifies the time scale of the BC forcing component of EBDragTempForcing relative to the time step size.

.. input_param:: EBDragTempForcing.use_idw_soil_temp_model

   **type:** Boolean, optional, default = false

   Enables the inverse distance weighting (IDW) soil temperature model, which computes local soil temperatures based on user-provided coordinates and temperatures.

.. input_param:: EBDragTempForcing.idw_x

   **type:** List of Reals, mandatory if IDW model is enabled

   The X-coordinates of the reference soil temperature points.

.. input_param:: EBDragTempForcing.idw_y

   **type:** List of Reals, mandatory if IDW model is enabled

   The Y-coordinates of the reference soil temperature points.

.. input_param:: EBDragTempForcing.idw_z

   **type:** List of Reals, mandatory if IDW model is enabled

   The Z-coordinates of the reference soil temperature points.

.. input_param:: EBDragTempForcing.idw_temp

   **type:** List of Reals, mandatory if IDW model is enabled

   The reference soil temperature values corresponding to each coordinate set.

.. input_param:: EBDragTempForcing.use_heatflux_model

   **type:** Boolean, optional, default = false

   Enables a heat flux-based soil temperature model where soil temperature is computed from surface temperature, heat flux, conductivity, and thickness.

.. input_param:: EBDragTempForcing.surface_heatflux

   **type:** Real, optional, default = 0.0

   Specifies the surface heat flux in W/m\ :sup:`2` \.

.. input_param:: EBDragTempForcing.thermal_conductivity

   **type:** Real, optional, default = 0.026

   Specifies the thermal conductivity of the soil/boundary layer in W/(m·K).

.. input_param:: EBDragTempForcing.boundary_layer_height

   **type:** Real, optional, default = 1.0

   Specifies the boundary layer height/thickness in meters for the heat flux model.