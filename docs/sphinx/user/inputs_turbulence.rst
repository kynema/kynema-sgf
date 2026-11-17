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
   

   

.. input_param:: OneEqKsgsM84_coeffs.density_weighted_production

   **type:** Boolean, optional, default = false

   Selects how the shear and buoyancy production terms of the subgrid TKE
   equation are built in the `OneEqKsgsM84` model. With the default
   (``false``) both terms are built per unit mass from
   :math:`\nu_t = \mu_t / \rho`, in the same way as the dissipation term,
   before the TKE equation multiplies the whole source by density. With
   ``true`` the terms are built from :math:`\mu_t` and carry an extra factor
   of :math:`\rho` relative to dissipation. This was the behavior before
   this option was added. It is kept only to reproduce earlier results.
   The two settings give the same answer when ``incflo.density = 1``.

.. input_param:: KLAxell_coeffs.density_weighted_production

   **type:** Boolean, optional, default = false

   Same as :input_param:`OneEqKsgsM84_coeffs.density_weighted_production`,
   for the shear and buoyancy production terms of the `KLAxell` model.
