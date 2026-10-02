.. _inputs_enst:

Section: Enstrophy
~~~~~~~~~~~~~~~~~~

This section controls  enstrophy  post processing. 
The prefix is the label set in ``incflo.post_processing``. For example
``incflo.post_processing = enst``


.. input_param:: enst.type

   **type:** String, mandatory

   To use enstrophy output specify with keyword ``Enstrophy``
   
.. input_param:: enst.output_interval

   **type:** Integer, optional, default = 10

   Specify the output frequency (in timesteps) for integrating enstrophy and 
   writing to file.
