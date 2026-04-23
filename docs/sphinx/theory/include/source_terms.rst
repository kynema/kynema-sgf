.. _source_terms:

Source terms
------------

Gravity Forcing
~~~~~~~~~~~~~~~~

The implementation of this source term allows the user to choose the full gravity term (:math:`\rho g`) or a perturbational form (:math:`(\rho - \rho_0) g`). By default, the full term is used, but the perturbational form can be turned on by adding ``ICNS.use_perturb_pressure = true`` to the input file.

The reference density (:math:`\rho_0`) is defined as ``1.0`` by default, can be defined as a constant through the input argument, ``incflo.density``, or can be defined as a spatially varying field within the flow setup (see physics/multiphase/Multiphase.cpp).

Using the perturbational form implies that the hydrostatic pressure is removed from the pressure variable, including its output. This means that the solution to the Poisson equation is actually the perturbational pressure, :math:`p'`, not :math:`p`. If the full pressure, :math:`p`, is desired for analysis or postprocessing purposes, the hydrostatic pressure can be added back to the pressure field via the input argument ``ICNS.reconstruct_true_pressure = true``. In order for this to operate in the code, the reference pressure field must be defined for the specific flow case being run. 

- An example of this is in physics/multiphase/Multiphase.cpp. To construct the reference pressure field, the reference gravity term must be integrated. This particular example assumes that the reference density only varies in z (or is constant), gravity acts only in z, and the hydrostatic pressure at the high z boundary is equal to 0. 

- In mathematical form, the derivation and calculation of the full pressure is as follows:

.. math:: \nabla p = \nabla p' + \rho_0 \boldsymbol{g}

- assume :math:`\boldsymbol{g} = g\hat{k}` and :math:`\frac{dp_0}{dz} = g\hat{k}`

.. math:: p = p' + \int_{z_{min}}^z \rho_0 g dz + p(z = z_{min}) 

- change reference frame to the top boundary, and assume :math:`p(z = z_{max}) = 0`
   
.. math:: p = p' - \int_z^{z_{max}} \rho_0 g dz + p(z = z_{max}) = p' - \int_z^{z_{max}} \rho_0 g dz

.. _mesoscale_forcing:

Mesoscale Forcing
~~~~~~~~~~~~~~~~~

To incorporate larger-scale atmospheric dynamics under real conditions,
Kynema-SGF offers two approaches. If mesoscale momentum and/or temperature
source terms are known exactly, e.g., from a numerical weather prediction (NWP)
model, then these may be directly applied. These mesoscale source terms would
come from the RHS of the mesoscale equations of motion and may also include the
effects of additional modeled physics such as radiation or moisture. This
mesoscale forcing approach is called the "tendencies" (or "mesoscale budget
components") approach. For more information, see `Draxl et al. (BLM 2021)
<https://doi.org/10.1007/s10546-020-00584-z>`_

If the mesoscale source terms are not known a priori, they may be derived on
the fly with a profile assimilation technique. This is an engineering approach
that applies a proportional controller to drive the instantaneous planar
averaged wind and/or temperature profiles towards known time--height data. This
approach can be used with NWP model output or observational data. For more
information, see `Allaerts et al. (BLM 2020)
<https://doi.org/10.1007/s10546-020-00538-5>`_

The application of these forcing approaches is detailed :ref:`here <inputs_meso_forcing>`.

Actuator Forcing
~~~~~~~~~~~~~~~~

Calculating actuator forces relies on sampling the velocity field at actuator points
at the beginning of each time step (*n*). Actuator-based models, i.e., actuator lines
and actuator disks, rely on internal implementations (e.g., Joukowsky disk, actuator-line wing)
or external turbine tools (OpenFAST) that use these sampled velocities to calculate forces 
and the motion of actuator points.
When the Godunov method is used, the motion of actuator points must be incorporated
into the application of actuator forces. This is because the Godunov method discretizes source terms
at the half time step (*n+1/2*). Therefore, the actuator force vectors are calculated using 
fluid velocities at *n*, and these actuator forces are applied at locations corresponding
to *n+1/2*.

