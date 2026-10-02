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

.. input_param:: DragTempForcing.blank_follow_fluid

   **type:** Boolean, optional, default = false (true with
   :input_param:`TerrainDrag.wall_treatment` ``= improved``)

   Relaxes the blanked (terrain) cells toward the temperature of the cell
   above them instead of ``DragTempForcing.soil_temperature``.
   Held at the soil temperature under air of another temperature, the
   blanked column conducts heat through its top face into the drag cell, a
   surface heat flux on top of the one the drag-cell forcing sets. With this
   option no gradient forms across the blanked face. The relaxation is
   integrated exactly over the time step,

   .. math::

      S = -\frac{1 - e^{-C_d \Delta t}}{\Delta t}\,(\theta_k - \theta_{k+1}),
      \qquad
      C_d = \min\left(\frac{c_d}{\Delta z\,|\mathbf{u}|}, \frac{10}{\Delta z}\right),

   with :math:`c_d` = :input_param:`DragTempForcing.drag_coefficient`; the
   explicit rate reaches 2 per step at :math:`\Delta z` = 5 m and
   :math:`\Delta t` = 1 s. The drag cells are not changed. The target is the
   cell directly above, so on slopes a blanked cell next to the air on its
   side still follows the cell above it.


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