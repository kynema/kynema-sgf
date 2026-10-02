.. _inputs_ke:
  
Section: KineticEnergy
~~~~~~~~~~~~~~~~~~~~~~

This section controls  kinetic energy  post processing. 
The prefix is the label set in ``incflo.post_processing``. For example
``incflo.post_processing = ke``


.. input_param:: ke.type

   **type:** String, mandatory

    To use kinetic energy post processing specify with keyword ``KineticEnergy``

.. input_param:: ke.output_interval

   **type:** Integer, optional, default = 10

   Specify the output frequency (in timesteps) for integrating kinetic energy
   and writing to file
