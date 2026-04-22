.. _inputs_ocean_waves:

Section: Ocean Waves
~~~~~~~~~~~~~~~~~~~~

This section is for setting up wave forcing and relaxation zones.

.. input_param:: OceanWaves.label

   **type:** String

   Name for Ocean Waves instance. As a placeholder, the label 'label' will be used in this documentation.

.. input_param:: OceanWaves.label.type

   **type:** String

   Type of wave to be used. Options include LinearWaves, StokesWaves, and W2AWaves. 
   Linear and stokes waves are monochromatic waves. W2AWaves uses the Waves2AMR library to read 
   and transform HOS-Ocean wave data into the Kynema-SGF domain, and this option requires the Waves2AMR 
   library (abbreviated W2A) to be enabled at compile time.

.. input_param:: OceanWaves.label.relax_zone_gen_length

   **type:** Real, optional, default = 4.0

   Length of region, in x-direction, over which the wave profile will be forced. 
   The region begins at the low x boundary and extends in the positive x direction, covering
   the entire width of the domain in the y direction and applying to all liquid at all depths
   (in z). Always on, even if argument is not specified.

.. input_param:: OceanWaves.label.numerical_beach_length

   **type:** Real, optional, default = 8.0

   Length of region, in x-direction, over which the waves will be forced to a flat interface. 
   The region begins at the high x boundary and extends in the negative x direction, covering
   the entire width of the domain in the y direction and applying to all liquid at all depths
   (in z). Always on, unless relax_zone_out_length is included in the input file.

.. input_param:: OceanWaves.label.numerical_beach_length_factor

   **type:** Real, optional, default = 1.0

   Length factor for numerical beach region. Following the formulation suggested by `Chen, Kelly, and Zang 
   (2019) <https://www.sciencedirect.com/science/article/pii/S0029801819303919>`_, this factor increases the 
   effective length of the numerical beach while keeping it as the same physical length. This approach stems
   from the observation that waves are typically fully absorbed in the first half of the numerical beach,
   leaving the second portion basically unused. Therefore, using a numerical beach length factor can
   reduce the size of the required domain while still achieving the same results. For example, if a numerical beach
   length of 2 wave lengths is needed, ``numerical_beach_length`` could be set to 1 wave length, and then
   ``numerical_beach_length_factor`` would be set to 2.

.. input_param:: OceanWaves.label.relax_zone_out_length

   **type:** Real, optional, default = 8.0

   Length of region, in x-direction, over which the wave profile will be forced. This input argument
   replaces, and thereby turns off, the numerical beach.
   The region begins at the high x boundary and extends in the negative x direction, covering
   the entire width of the domain in the y direction and applying to all liquid at all depths
   (in z). Only on when included in the input file.

.. input_param:: OceanWaves.label.zero_sea_level

   **type:** Real, optional, default = 0.0

   Default location (in z) of the liquid surface. Where the water would rest in the absence of waves.

.. input_param:: OceanWaves.label.water_depth

   **type:** Real, optional, default = 0.5

   The depth of the water in the simulation, used for calculating analytical wave profiles (linear 
   and stokes waves). Should be equal to the difference between the zero_sea_level and the low z boundary location.

.. input_param:: OceanWaves.label.timeramp_period

   **type:** Real, optional

   An initial ramp-up period for the wave forcing. Without specifying a period, the wave 
   forcing will begin at full strength.

.. input_param:: OceanWaves.label.initialize_wave_field

   **type:** Boolean, optional, default = false

   By default, the domain will be initialized with a flat interface; if this option
   is turned on, the wave profile will be initialized over the entire domain. If there is a specified
   relax_zone_out_length, this option is automatically turned on.

.. input_param:: OceanWaves.label.current

   **type:** Real, optional, default = 0.0

   Flow speed of a uniform current in the x-direction to be added to the wave forcing. This option is
   currently only compatible with the LinearWaves wave type and the numerical beach. If a nonzero current is
   specified for a case with nonlinear waves, the code will abort.

The following input arguments are only valid for the LinearWaves and StokesWave wave types:

.. input_param:: OceanWaves.label.wave_length

   **type:** Real, mandatory

   The wave length of the wave profile. This argument can be omitted for Stokes waves if
   the wave period is provided instead.

.. input_param:: OceanWaves.label.wave_height

   **type:** Real, mandatory

   The amplitude of the wave profile

.. input_param:: OceanWaves.label.wave_phase_offset_radians

   **type** Real, optional, default = 0.

   This specifies a phase offset to the target wave profile.
   It is in terms of phase (radians) relative to the wave, not
   in terms of time. Note that this argument should be retained
   when restarting a simulation so that the phase of the wave forcing
   is consistent with the preceding simulation.

.. input_param:: OceanWaves.label.wave_phase_offset_degrees

   **type** Real, optional, default = 0.

   This specifies a phase offset to the target wave profile
   in terms of degrees. If the offset is specified in both radians
   and degrees, the code will abort; only one of these arguments
   can be used at a time.

The following input arguments are only valid for the StokesWave wave type:

.. input_param:: OceanWaves.label.order

   **type:** Integer, mandatory

   The order of the Stokes wave formula being used. All Stokes wave theory (wave profile and
   dispersion relation) is applied from `Fenton (1985)
   <https://ascelibrary.org/doi/10.1061/%28ASCE%290733-950X%281985%29111%3A2%28216%29>`_.
   The minimum order is 2, and the maximum order is 5.

.. input_param:: OceanWaves.label.wave_period

   **type:** Real, optional

   If the wave period is provided and the wave length is not, the wave length will be solved
   iteratively using the wave height, the wave period, and the Stokes waves dispersion relation.
   If the wave length is not provided, this argument becomes mandatory.

.. input_param:: OceanWaves.label.stokes_wavelength_order

   **type:** Real, optional

   Specifies the order of the dispersion relation used to calculate the wave length. By default,
   this is equal to the ``order`` of the waves. Practically, the minimum value is 1, and the
   maximum is 5.

.. input_param:: OceanWaves.label.stokes_wavelength_tolerance

   **type:** Real, optional, default = 1e-10

   Convergence tolerance of the iterative process to calculate the wave length.

.. input_param:: OceanWaves.label.stokes_wavelength_iter_max

   **type:** Integer, optional, default = 40

   Maximum number of iterations during the process to calculate the wave length.

The following input arguments are only valid for the W2AWaves wave type:

.. input_param:: OceanWaves.label.HOS_modes_filename

   **type:** String, mandatory

   The name of the modes file, output by HOS-Ocean or HOS-NWT, in the SWENSE format.

.. input_param:: OceanWaves.label.HOS_simulation_is_ocean

   **type:** Boolean, optional, default = true

   Specify whether the modes file is from HOS-Ocean or HOS-NWT. If this is specified incorrectly,
   the code will abort by detecting mismatched parameters in the modes file.

.. input_param:: OceanWaves.label.HOS_init_timestep

   **type:** Integer, optional, default = 0

   The time step in the modes file for the Kynema-SGF simulation to start at.

.. input_param:: OceanWaves.label.HOS_init_time

   **type:** Integer, optional, default = 0

   The physical time in the modes file for the Kynema-SGF simulation to start at.
   This argument is only active if HOS_init_timestep is omitted. Kynema-SGF will pick the
   time step in the modes closest to the specified time.

.. input_param:: OceanWaves.label.HOS_domain_offset_x

   **type:** Real, optional, default = 0.

   A physical space offset for mapping between the HOS domain and the Kynema-SGF domain.
   The coordinates of the HOS domain always start at (x,y) = (0,0), but often the setup
   of an Kynema-SGF domain is not constrained in this way, and also the Kynema-SGF domain
   may only represent a portion of the HOS domain. This parameter, :math:`x_textrm{off}`, 
   contributes to the mapping such that :math:`x_\textrm{Kynema-SGF} = x_\textrm{HOS} + x_\textrm{off}`.

.. input_param:: OceanWaves.label.HOS_domain_offset_y

   **type:** Real, optional, default = 0.

   A physical space offset for mapping between the HOS domain and the Kynema-SGF domain.
   This parameter, :math:`y_textrm{off}`, contributes to the mapping such that
   :math:`y_\textrm{Kynema-SGF} = y_\textrm{HOS} + y_\textrm{off}`.

.. input_param:: OceanWaves.label.fftw_planner_flag

   **type:** String, optional, default = estimate

   When setting up a plan for the inverse Fourier transform within the Waves2AMR library,
   the FFTW algorithm can use different techniques to choose among available methods. Some of these
   are faster than others, and the optimal choice can also depend on the architecture. The
   default, "estimate", which corresponds to FFTW_ESTIMATE, is deterministic. The other
   options are "exhaustive", "patient", and "measure". Variations from nondeterministic
   approaches are tiny, on the order of machine precision.

.. input_param:: OceanWaves.label.number_interp_points_in_z

   **type:** Integer, mandatory

   When Waves2AMR converts mode data to spatial data, the z locations to perform the transformation
   must be chosen. Kynema-SGF does this by creating a geometric series from the water surface
   to the lower z boundary, starting with a small spacing between points near the surface and 
   then expanding downward. This argument is the total number of z locations where wave modes will be
   transformed to inform the wave velocity field. After the transformation step, the data is then interpolated
   to the Kynema-SGF mesh. More points means better resolution over the whole depth, but more points
   also means more work for the solver.

.. input_param:: OceanWaves.label.interp_spacing_at_surface

   **type:** Real, mandatory

   The physical spacing between interpolation points at the water surface.
   This is the most influential parameter for the geometric series that dictates the location of the points in z.
   This should be set to near the mesh spacing in z around the water surface.

.. input_param:: OceanWaves.label.number_interp_above_surface

   **type:** Integer, optional, default = 1

   The number of points placed above the mean water surface for the velocity transformation process. The spacing
   between the points above the surface is equal to the ``interp_spacing_at_surface``. When setting this value, the wave height
   should be considered so that velocity can be accurately computed for portions of the waves above the mean surface.
