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

.. input_param:: OneEqKsgsM84.surfaceRANS

   **type:** Boolean, optional, default = false

   Blends the filter width of the `OneEqKsgsM84` model with a wall-distance
   length scale near the ground, in the same way as the `surfaceRANS` option
   of the `Kosovic` model. The length scale of a refined box that touches the
   ground then no longer depends on the refinement level, which removes the
   slow near-wall layer and the velocity bands along the box faces (issue
   1806). The blended length scale is used for the eddy viscosity, the
   turbulent diffusivity and the dissipation of the subgrid kinetic energy:
   :math:`l^2 = (1-f)^n \Delta^2 + f^n l_{RANS}^2` with
   :math:`f = \exp(-z / z_s)` and
   :math:`l_{RANS} = \kappa z / \phi_m \,(C_\epsilon / C_e^3)^{1/4}`,
   the length for which the model recovers the log-law eddy viscosity
   :math:`\kappa u_* z`. The stability function :math:`\phi_m` uses
   ``ABL.monin_obukhov_length`` when it is given and otherwise the Obukhov
   length of the ABL wall function at the current step.

.. input_param:: OneEqKsgsM84.switchLoc

   **type:** Real, optional, default = 24.0

   Height :math:`z_s` (m) over which the blending weight decays when
   `OneEqKsgsM84.surfaceRANS` is active. About four to five vertical cells of
   the level surrounding the refined box keeps the length scale at or above
   the coarse value through the layer where the box cannot resolve its own
   stress; the filter width is recovered beyond about five times this height.

.. input_param:: OneEqKsgsM84.surfaceRANSExp

   **type:** Real, optional, default = 2.0

   Exponent :math:`n` applied to the blending weights when
   `OneEqKsgsM84.surfaceRANS` is active. 
   

   
