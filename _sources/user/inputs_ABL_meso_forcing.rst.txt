.. spelling:word-list::

   ntime
   nheight
   tflux
   directToConstant
   indirectToConstant
   indirectToDirect

.. _inputs_meso_forcing:

Section: ABLMesoForcing
~~~~~~~~~~~~~~~~~~~~~~~

This section is for setting mesoscale forcing parameters for
mesoscale-to-microscale coupling (MMC). With the exception of
``ABL.mesoscale_forcing`` and ``ABL.tendency_forcing``, the parameters below
apply to both ``ABLMesoForcingMom`` and ``ABLMesoForcingTemp`` for momentum and
temperature, respectively.

.. input_param:: ABL.mesoscale_forcing

   **type:** String, optional

   Path to a netcdf file with mesoscale data. This datafile should have
   dimensions "ntime" and "nheight". Time-height-varying variables include
   "momentum_u", "momentum_v", and "temperature" (virtual potential
   temperature); time-varying variables include "tflux" (kinematic heat flux)
   and "transition_height" (optional, only used if
   :input_param:`ABLMesoForcing*.forcing_transition` is specified without a
   :input_param:`ABLMesoForcing*.constant_transition_height`). The time-height
   data are stored as flattened 1-D arrays, wherein the height dimension
   changes the fastest.

.. input_param:: ABL.tendency_forcing

   **type:** Boolean, optional, default = false
   
   If "true", momentum and temperature source terms are directly provided
   within the :input_param:`ABL.mesoscale_forcing` netcdf file; 
   :input_param:`ABLMesoForcing*.forcing_scheme` is ignored. Typically, these
   source terms are determined from mesoscale model tendencies. 
   
   If "false", the input velocity and temperature time-height data are used for
   profile assimilation forcing schemes.
   
.. input_param:: ABLMesoForcing*.forcing_scheme

   **type:** String, mandatory

   Options include: "direct" or "indirect" profile assimilation. **Direct
   profile assimilation (DPA)** is a strong enforcement of the input mesoscale
   profiles. Height-varying source terms are calculated at every time step to
   drive the error between the instantaneous planar-averaged fields and the
   mesoscale data to zero. This strategy implies a high level of confidence in
   the specified mesoscale conditions and can result in excessive turbulence
   resolved in the microscale.

   In comparison, **indirect profile assimilation** weakly enforces the
   mesoscale profiles by using an approximation to the DPA forcing. This is an
   engineering approach that allows for local deviation from the mesoscale mean
   quantities and allows the microscale flow to find its own equilibrium
   between mean profiles and the associated resolved turbulence. The
   approximate forcing is derived by applying a polynomial regression to the
   DPA forcing profile.

.. input_param:: ABLMesoForcing*.control_gain

   **type:** Real, optional, default = 0.2

   A relaxation factor used to scale the magnitude of the source term, which
   may be interpreted as an inverse time scale with units of 1/s.

.. input_param:: ABLMesoForcing*.forcing_transition

   **type:** String, optional, default = "none"

   Used to specify two forcing schemes for upper and lower regions under more
   complicated MMC modeling scenarios, e.g., with complex background conditions
   and/or uncertainty in the mesoscale data. The default is apply the
   :input_param:`ABLMesoForcing*.forcing_scheme` for the full mesoscale profile
   over the entire microscale domain. Other options include "directToConstant",
   "indirectToConstant", and "indirectToDirect".
   
   When blending to a constant forcing profile, the slope of the forcing
   profile at the transition height (specified by
   :input_param:`ABLMesoForcing*.constant_transition_height` or provided in the
   :input_param:`ABL.mesoscale_forcing` datafile) is linearly reduced to 0 from
   the transition height up to the transition height +
   :input_param:`ABLMesoForcing*.transition_thickness`.

   When blending to the DPA forcing profile (e.g., in the free atmosphere), the
   two forcing profiles are linearly blended from one to the other over the
   :input_param:`ABLMesoForcing*.transition_thickness`, starting from the
   transition height as described above.


Indirect Profile Assimilation
-----------------------------
The following parameters are specific to the IPA scheme
(:input_param:`ABLMesoForcing*.forcing_scheme` = "indirect"). At the moment,
only third-order polynomial regression is supported.

.. input_param:: ABLMesoForcing*.weighting_heights

   **type:** List of Reals (has to be same length as
   :input_param:`ABLMesoForcing*.weighting_values`), optional

   Height(s) in meters at which IPA regression weights are provided.
   
.. input_param:: ABLMesoForcing*.weighting_values

   **type:** List of Reals (has to be same length as
   :input_param:`ABLMesoForcing*.weighting_heights`), optional

   IPA regression weights at the corresponding
   :input_param:`ABLMesoForcing*.weighting_heights`. The default behavior is to
   use uniform weighting. Nonuniform weighting is generally ill-advised as
   runaway positive or negative forcing values may be possible.

.. input_param:: ABLMesoForcing*.normalize_by_zmax

   **type:** Boolean, optional, default = false

   If "true", the height coordinate is normalized by the domain height when
   performing the IPA regression. Provided for consistency with a legacy solver
   implementation to improve conditioning of the regression matrix but should
   *not* be needed.


Partial Profile Assimilation
-----------------------------
The following parameters are for "partial" profile assimilation, enabled by
:input_param:`ABLMesoForcing*.forcing_transition` being not set to "none". This
will only partially apply the instantaneous IPA forcing profiles over the
simulation domain. Above a specified transition layer, a secondary forcing
profiles may be applied.

.. input_param:: ABLMesoForcing*.transition_thickness

   **type:** Real

   The thickness of the layer over which the forcing scheme transitions from
   the lower scheme to the upper scheme. 

.. input_param:: ABLMesoForcing*.constant_transition_height

   **type:** Real

   The base of the transition layer, which is invariant for the duration of the
   entirety of the simulation. To specify a time-varying transition layer
   height that, e.g., tracks the evolution of the ABL height, omit this
   parameter and include the time-varying ``transition_height`` variable within
   the :input_param:`ABL.mesoscale_forcing` datafile.

