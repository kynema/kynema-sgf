.. _inputs_boundary_conditions:

Section: Boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section controls the boundary conditions. Only non-periodic BC's need to be defined here.

.. input_param:: xlo.type (or ylo.type, zlo.type, xhi.type, yhi.type, zhi.type)

   **type:** String, mandatory if non-periodic

   Boundary condition type on the lo (or hi) side of the domain.
   Current options are: periodic, pressure_outflow, mass_inflow,
   mass_inflow_outflow, no_slip_wall, slip_wall, symmetric_wall and wall_model.

.. input_param:: xlo.temperature (or ylo.temperature, etc)

   **type:** Real, optional, default = 0.0

   Specifies a temperature gradient at the wall, only activated for slip_wall and wall_model BC types.

Mass inflow boundary conditions
```````````````````````````````

For the mass_inflow boundary condition type, the user can specify constant values in the input file::

  xlo.type = "mass_inflow"
  xlo.velocity = 1.0 1.0 0.0
  xlo.temperature = 300.0

In addition to being able to specify constant values, the user has access to a variety of included complex boundary conditions (e.g., LinearProfile, PowerLawProfile, Rankine) that are activated in the input file. For example::

  xlo.type = "mass_inflow"
  xlo.density = 1.0
  xlo.velocity = 5.0 5.0 0.0
  xlo.velocity.inflow_type = LinearProfile

LinearProfile
"""""""""""""

The linear profile uses the following input parameters. The linear profile can be use for multiple fields and therefore the input file option should contain the ``field_name``, e.g. ``LinearProfile.velocity.direction = 1``. It sets up the following profile.

.. math::

   f(x_i) =  \begin{cases}
      f_{\text{min}} & \text{if } x_i < x_{\text{min}}\\
      f_{\text{max}} & \text{if } x_i > x_{\text{max}}\\
      f_{\text{min}} + (f_{\text{max}}-f_{\text{min}}) \frac{x_i-x_{\text{min}}}{x_{\text{max}}-x_{\text{min}}} & \text{otherwise}
    \end{cases}

.. input_param:: LinearProfile.<field_name>.direction

   **type:** int, optional, default = 2

   The coordinate direction along which this profile should be initialized (wall normal direction). :math:`i` value in the equation.

.. input_param:: LinearProfile.<field_name>.start

   **type:** Real, optional, default = minimum of domain in the ``direction`` direction

   Start location of the linear profile. :math:`x_{\text{min}}` value in the equation.

.. input_param:: LinearProfile.<field_name>.stop

   **type:** Real, optional, default = maximum of domain in the ``direction`` direction

   Stop location of the linear profile. :math:`x_{\text{max}}` value in the equation.

.. input_param:: LinearProfile.<field_name>.start_val

   **type:** Real, required (length = number of components of ``field_name``)

   Start value of field. :math:`f_{\text{min}}` value in the equation.

.. input_param:: LinearProfile.<field_name>.stop_val

   **type:** Real, required (length = number of components of ``field_name``)

   Stop value of field. :math:`f_{\text{max}}` value in the equation.

Rankine
"""""""

The Rankine boundary condition uses the following input parameters. It is used for velocity fields only and the reference velocity is specified by ``incflo.velocity``.

.. input_param:: Rankine.Umax

   **type:** Real, optional, default = 15.0

   The maximum velocity.

.. input_param:: Rankine.Rmax

   **type:** Real, optional, default = 125.0

   The maximum radius.

.. input_param:: Rankine.start_location

   **type:** List of Real, optional, default = [-10 * ``Rmax``, 0.0, 0.0]

   The start location.


PowerLawProfile
"""""""""""""""

The power law profile uses the following input parameters. This one is only for velocity so there is no need to insert the field name in the input line. This is the equation describing this boundary condition.

.. math::

   \begin{align}
   p_f &= \text{min}\left(\text{max}\left(\left(\frac{x_i - x_o}{x_r}\right)^s, \frac{u_{\text{min}}}{u_r}\right), \frac{u_{\text{max}}}{u_r}\right) \\
   f(x_i, n) &= u_n  p_f
   \end{align}

.. input_param:: PowerLawProfile.direction

   **type:** int, optional, default = 2

   The coordinate direction along which this profile should be initialized (wall normal direction). :math:`i` value in the equation.

.. input_param:: PowerLawProfile.zref

   **type:** Real, required

   The reference height. :math:`x_r` in the equation.

.. input_param:: PowerLawProfile.shear_exponent

   **type:** Real, required

   The shear exponent value. :math:`s` in the equation.

.. input_param:: PowerLawProfile.uref

   **type:** Real list, required

   The reference value of the velocity vector used to propagate the plane. :math:`u_r` in the equation.

.. input_param:: PowerLawProfile.zoffset

   **type:** Real, optional, default = 0.0

   The offset in the direction (should be zero for a boundary condition). :math:`x_o` in the equation.

.. input_param:: PowerLawProfile.umin

   **type:** Real, default = 0.0

   The minimum velocity cutoff in the mean power law profile. :math:`u_{\text{min}}` in the equation.

.. input_param:: PowerLawProfile.umax

   **type:** Real, default = 0.0

   The maximum velocity cutoff in the mean power law profile. :math:`u_{\text{max}}` in the equation.

Custom boundary conditions
""""""""""""""""""""""""""

If the user wants to define their own boundary conditions, this is done by editing `CustomScalar` and `CustomVelocity` source and header files in the `udfs` folder. `CustomScalar` is used for scalar fields and `CustomVelocity` is used for velocity fields. These can then be activated in the input file as such::

  xlo.type = "mass_inflow"
  xlo.temperature.inflow_type = CustomScalar
  CustomScalar.foo = 1.0
  xlo.velocity.inflow_type = CustomVelocity
  CustomVelocity.foo = 1.0

They do not both need to be defined at the same time. It is the user's responsibility to ensure that the source files are appropriately edited for their use case. Examples of how these files can be edited are found through comparison of the other mass_inflow functions in the `udfs` folder.

Mass inflow-outflow boundary conditions
```````````````````````````````````````

The mass_inflow_outflow boundary condition is designed to handle both inflow and outflow at the same boundary.
For the advection schemes, it implements a Neumann type behavior at the outflow cells and a Dirichlet behavior at the inflow cells.
It uses Neumann conditions for the MAC and nodal projections and
enforces solvability before the projections
by correcting the outflow to match with the inflow within the specified mass_inflow_outflow boundaries.
It uses a Dirichlet condition for the diffusion solver.

Both the approaches mentioned above for the mass inflow condition,
constant values and UDFs, can be used to specify the boundary values.
The outflow values will be automatically replaced by a value from the interior cell
to enforce the Neumann type behavior.
See the ``freestream_godunov_inout`` test for an example that uses the TwoLayer UDF.
This test involves two z-layers of the flow along opposite x-directions.
The input file options are copied here::

  geometry.is_periodic  =  0   1   0   # Periodic in y

  # Boundary conditions
  TwoLayer.bottom_vel   = -1.0 0.0 0.0
  TwoLayer.top_vel      =  1.0 0.0 0.0
  TwoLayer.init_perturb = 0.9
  TwoLayer.z_partition  = 0.5

  xlo.type = "mass_inflow_outflow"
  xlo.density = 1.0
  xlo.velocity.inflow_outflow_type = TwoLayer

  xhi.type = "mass_inflow_outflow"
  xhi.density = 1.0
  xhi.velocity.inflow_outflow_type = TwoLayer

  zlo.type = "slip_wall"
  zhi.type = "slip_wall"


The most applicable use case for this boundary condition is with the
:ref:`kynemasgf-abl-bndry-io` for flows that change directions
across the vertical coordinate or with time.
See the ``abl_bndry_input_native_inout`` test for an example,
the relevant inputs are included below.

.. code-block:: none

  ABL.bndry_file = "../abl_bndry_output_native/bndry_files"
  ABL.bndry_io_mode = 1
  ABL.bndry_var_names = velocity temperature
  ABL.bndry_output_format = native

  xlo.type = "mass_inflow_outflow"
  xlo.density = 1.0
  xlo.temperature = 0.0   # dummy

  xhi.type = "mass_inflow_outflow"
  xhi.density = 1.0
  xhi.temperature = 0.0

  ylo.type = "mass_inflow_outflow"
  ylo.density = 1.0
  ylo.temperature = 0.0

  yhi.type = "mass_inflow_outflow"
  yhi.density = 1.0
  yhi.temperature = 0.0

If the inflow files do not include a specific field or plane,
then the relevant BC can be specified in the input file.
Note that if unspecified, the velocities at a boundary would default to zero,
which implements a Dirichlet type behavior.
If an outflow is expected at that boundary, any non-zero outflow velocity may be specified,
which would implement a Neumann type behavior as the
boundary cell velocity will be overwritten by the inner cell velocity.

TwoLayer
""""""""

The two layer boundary condition uses the following input options. It is only used for velocity.

.. input_param:: TwoLayer.top_vel

   **type:** List of Real, required

   The velocity in the top of the domain.

.. input_param:: TwoLayer.bottom_vel

   **type:** List of Real, required

   The velocity in the bottom of the domain.

.. input_param:: TwoLayer.init_perturb

   **type:** Real, optional, default = 1.0

   The initial perturbation.

.. input_param:: TwoLayer.z_partition

   **type:** Real, optional, default = 0.5

   The z-coordinate for the top and bottom.


Dynamic wall model (Wave model)
```````````````````````````````
The Moving Surface Drag (MOSD) model developed by `Ayala et al (2024) <https://doi.org/10.1007/s10546-024-00884-8>`_ is used as the dynamic wall model. The model calculates the stress (form drag) imparted by a moving wave. The model enables wave phase-resolving physics without the use of wave-phase adapting computational grids.

.. input_param:: wave_mosd.amplitude
   **type:** Real, required, default = 0.05

   Specifies the amplitude of the wave, only activated if ``WallFunction.wall_shear_stress_type = mosd``

.. input_param:: wave_mosd.wavenumber
   **type:** Real, required, default = 4

   Specifies the wave number of the wave, only activated if ``WallFunction.wall_shear_stress_type = mosd``

.. input_param:: wave_mosd.frequency
   **type:** Real, required, default = 0.8

   Specifies the frequency of the wave, only activated if ``WallFunction.wall_shear_stress_type = mosd``

Example::

  zlo.type =   "wall_model"
  WallFunction.wall_shear_stress_type = mosd
  wave_mosd.amplitude = 0.05
  wave_mosd.wavenumber = 4
  wave_mosd.frequency = 0.8

.. note:: This wall model is only applicable for the lower boundary ``zlo.type``. Also, it is set for only monochromatic waves.

Currently, the dynamic wall model is only available for ``incflo.physics = ChannelFlow``. The work to integrate this condition with the ABL class is under progress. See the ``channel_mosd`` test for an example that uses the dynamic wall model.

Pressure outflow boundary conditions
````````````````````````````````````

The pressure_outflow boundary condition is the most common boundary condition used
when flow out of a boundary is desired. By default, this sets the pressure at the outflow
plane to 0 and assumes a zero gradient for other flow quantities (e.g., velocity and
scalars). Also by default, this boundary condition clips fluxes that would be advected
into the domain, which is for the sake of stability.

Having a uniform pressure value at the outflow is not often physically valid for
flows of interest, such as stratified ABLs and ocean waves. However, instead of changing
the target pressure of the boundary condition, it is typically more useful to change the
formulation of the source terms, transforming the pressure variable into the difference
between the true pressure and some constant, non-uniform reference pressure profile.
In the context of stratified ABLs, which typically apply gravity through the
BoussinesqBuoyancy source term, this pressure modification is realized with the
additional source term ABLMeanBoussinesq. In the context of ocean waves, which typically
apply gravity through the GravityForcing source term, this pressure modification is
realized through the input option "ICNS.use_perturb_pressure".

Finally, the default clipping of inflow at pressure_outflow boundaries can be disabled.
This is not recommended, but it is possible with the input option
"allow_inflow_at_pressure_outflow". This input argument is appended to the PDE name
where it should be applied (e.g., ICNS, temperature, or tke).
