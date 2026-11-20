.. _inputs_turbulence:

Section: turbulence
~~~~~~~~~~~~~~~~~~~

This section is for setting turbulence model parameters

.. input_param:: turbulence.model

   **type:** String, optional, default = Laminar

   Specifies which turbulence model to use, by default "Laminar" is
   chosen (effectively no turbulence model).  Currently the supported
   turbulence models are "Smagorinsky", "AMD", "Kosovic", 
   "OneEqKsgsM84", "KOmegaSST", "KOmegaSSTIDDES" or "KLAxell".

   
.. input_param:: Smagorinsky_coeffs.Cs

   **type:** Real, optional, default = 0.135

   Specifies the coefficient used in the `Smagorinsky` turbulence model. 
   

   
