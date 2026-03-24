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