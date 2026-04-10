.. _inputs_abl:

Section: ABL
~~~~~~~~~~~~

This section is for setting atmospheric boundary layer parameters.

.. input_param:: ABL.kappa

   **type:** Real, optional, default = 0.41

   Wall model coefficient.

.. input_param:: ABL.surface_roughness_z0

   **type:** Real, optional, default = 0.1

   Wall model surface roughness length (in meters). When specified, this sets both the aerodynamic 
   (pertaining to velocity) and thermal (pertaining to temperature) roughness lengths.
   These roughness lengths can be set separately using the two parameters listed below.
   Specifying these separately will override the combine parameter; however,
   if this parameter (which sets both) is used with only one of the roughness parameters below,
   the code will abort due to the ambiguity.

.. input_param:: ABL.aerodynamic_roughness_length

   **type:** Real, optional, default = 0.1

   Wall model surface roughness length (in meters) for aerodynamic processes, i.e., the wall momentum flux/shear stress model.

.. input_param:: ABL.thermal_roughness_length

   **type:** Real, optional, default = 0.1

   Wall model surface roughness length (in meters) for thermal processes, i.e., the wall temperature flux model.

.. input_param:: ABL.normal_direction

   **type:** Integer, optional, default = 2

   Wall model normal direction. x-direction = 0, y-direction = 1, z-direction = 2.

.. input_param:: ABL.log_law_height

   **type:** Real, optional

   Height to evaluate the log law for the wall model.
   Currently, if this parameter is not specified in the input file, the first half cell height is calculated
   and is used to set the log law height.
   Therefore the log law height depends on the domain size and number of elements in the normal direction.
   If this parameter is set to a number the log law is evaluated at a fixed height.
   Note: currently the fluctuating velocity terms in the shear stress model are only
   available at the first cell center above the wall. This limitation will be removed soon.

.. input_param:: ABL.surface_temp_flux

   **type:** Real, optional

   Surface temperature flux setting for the ABL wall function. Specifies a constant temperature flux
   at the wall-modeled boundary. This is not a required argument because there are other options for
   setting up the surface temperature condition.

.. input_param:: ABL.surface_temp_timetable

   **type:** String, optional

   File name of surface temperature time table, allowing the surface temperature
   to change with time without specifying a surface temperature rate.

.. input_param:: ABL.surface_temp_rate

   **type:** Real, optional

   Constant rate at which the surface temperature changes.

.. input_param:: ABL.surface_temp_init

   **type:** Real, optional

   Initial temperature of the wall-modeled surface. This parameter is only active
   when the surface temperature rate is specified. If this parameter is active but not
   specified, the initial temperature will be set to the reference temperature of the simulation.

.. input_param:: ABL.surface_temp_rate_tstart

   **type:** Real, optional

   Start time of the surface temperature rate of change. Prior to this time,
   the surface temperature remains at the initial value. This parameter is only
   active when the surface temperature rate is specified. The default start time is 0.


.. input_param:: ABL.temperature_heights

   **type:** List of Reals, mandatory

   Height(s) in meters at which temperature values are prescribed.

.. input_param:: ABL.temperature_values

   **type:** List of Reals (has to be same length as :input_param:`ABL.temperature_heights`), mandatory

   Temperature values in Kelvin at the corresponding :input_param:`ABL.temperature_heights`.
   The temperature below the first height is assumed to be constant and equal to the
   first temperature value.
   The temperature between values is initialized to have linear variation.
   The final temperature is constant above the last specified height.


.. input_param:: ABL.perturb_velocity

   **type:** Boolean, optional, default = true

   If true this flag turns on perturbations to the freestream flow.

.. input_param:: ABL.pertub_ref_height

   **type:** Real, optional, default = 50.0

   Reference height for velocity perturbations,
   perturbations exist below this height and decay above this height.
   Only active when :input_param:`ABL.perturb_velocity` = true.

.. input_param:: ABL.Uperiods

   **type:** Real, optional, default = 4.0

   Number of sinusoidal waves in x-direction.
   Only active when :input_param:`ABL.perturb_velocity` = true.

.. input_param:: ABL.Vperiods

   **type:** Real, optional, default = 4.0

   Number of sinusoidal waves in y-direction.
   Only active when :input_param:`ABL.perturb_velocity` = true.

.. input_param:: ABL.deltaU

   **type:** Real, optional, default = 1.0

   Amplitude of fluctuations in x-direction.
   Only active when :input_param:`ABL.perturb_velocity` = true.

.. input_param:: ABL.deltaV

   **type:** Real, optional, default = 1.0

   Amplitude of fluctuations in y-direction.
   Only active when :input_param:`ABL.perturb_velocity` = true.

.. input_param:: ABL.perturb_temperature

   **type:** Boolean, optional, default = false

   Perturb temperature field with random fluctuations.

.. input_param:: ABL.theta_amplitude

   **type:** Real, optional, default = 0.8 K

   Amplitude of the temperature perturbations added to the initial field. Only
   active when :input_param:`ABL.perturb_temperature` is true.

.. input_param:: ABL.cutoff_height

   **type:** Real, optional, default = domain height

   Height below which temperature perturbations are added

.. input_param:: ABL.random_gauss_mean

   **type:** Real, optional, default = 0.0

   Mean for the Gaussian random number generator

.. input_param:: ABL.random_gauss_var

   **type:** Real, optional, default = 1.0

   Variance for the Gaussian random number generator


.. input_param:: ABL.bndry_file

   **type:** String, optional, default = ""

   NetCDF-4 file name for ABL inflow

.. input_param:: ABL.bndry_io_mode

   **type:** Int, optional, default = -1

   IO mode (0=output, 1=input)

.. input_param:: ABL.bndry_planes

   **type:** String, optional, default = ""

   IO planes for ABL inflow

.. input_param:: ABL.bndry_output_start_time

   **type:** Real, optional, default = 0.0

   Time at which to start ABL inflow output

.. input_param:: ABL.bndry_var_names

   **type:** String, optional, default = ""

   Variables for IO for ABL inflow

.. input_param:: ABL.wall_shear_stress_type

   **type:** String, optional, default = "Moeng"

   Wall shear stress model: options include
   "constant", "local", "Schumann", and "Moeng"

.. input_param:: ABL.bndry_output_format

   **type:** String, optional, default = "native"

   Output of boundary plane files. Valid values are ``netcdf`` and ``native``.

.. input_param:: ABL.initial_condition_input_file

   **type:** String, optional, default= ""

   File that contains initial conditions for the
   velocity field in netcdf file format.
   This file is expected to have the same dimensions as the simulation.
   Values are passed directly from the file to the velocity field inside the code.
   Only spanwise velocity components are supported.

.. input_param:: ABL.anelastic

   **type:** Boolean, optional, default= false

   Activate anelastic behavior. This adds `reference_density` and
   `reference_pressure` fields.

.. input_param:: ABL.bottom_reference_pressure

   **type:** Real, optional, default = 1.01325e5

   Reference pressure at the bottom of the domain. Used for anelastic ABL.

.. input_param:: ABL.initial_wind_profile 

   **type:** Boolean, optional, default= false

   Activates the reading of wind speed profile from a file. Recommended for 
   RANS models and also for wind conditions input from climate model.

.. input_param:: ABL.rans_1dprofile_file 

   **type:** String, optional, default = ""

   This input is required when the ABL.initial_wind_profile is set to True. 

.. input_param:: ABL.minimum_vertical_position 

   **type:** Real, optional, default = 0.1

   When using an initial wind profile, this parameter limits how small the vertical
   position can be. It affects calculations in the initialization and in DragForcing.

.. input_param:: ABL.meso_sponge_start 

   **type:** Real, optional, default = 600

   Approximate height of the planetary boundary layer height to enable the forcing 
   in the free atmosphere. Recommended for use with RANS model and optionally to run 
   LES with non canonical flow conditions. If not specified in the input file, the 
   default value is active for the TKE source term but not for velocity or temperature.
   For this value to be used, the source terms must be included as input arguments 
   as follows, for temperature, velocity, or TKE, respectively.
   
   `Temperature.source_terms  = TemperatureFreeAtmosphereForcing`

   `ICNS.source_terms  = VelocityFreeAtmosphereForcing`

   `TKE.source_terms = KransAxell`

.. input_param:: ABL.velocity_sponge_start

   **type:** Real, optional

   Instead of specifying the same sponge start location for velocity, temperature, and TKE,
   this input argument allows the velocity sponge start (involved in `VelocityFreeAtmosphereForcing`)
   to be differentiated from the value specified by :input_param:`ABL.meso_sponge_start`.

.. input_param:: ABL.temp_sponge_start

   **type:** Real, optional

   Instead of specifying the same sponge start location for velocity, temperature, and TKE,
   this input argument allows the temperature sponge start (involved in `TemperatureFreeAtmosphereForcing`)
   to be differentiated from the value specified by :input_param:`ABL.meso_sponge_start`.

.. input_param:: ABL.wall_het_model

   **type:** String, optional, default = "none"

   Allows the use of different surface model options for the Monin-Obukhov length. Currently supports two options:
   (i) "none" - original model in the code and (ii) "mol" - Monin-Obukhov length is constant while heat-flux varies 

.. input_param:: ABL.monin_obukhov_length

   **type:** Real, optional, default = -1e30 

   Used in conjunction with `ABL.wall_het_model`. The default value runs a neutral boundary layer. 

.. input_param:: ABL.terrain_aligned_profile 

   **type:** Boolean, optional, default= false

   Used in conjunction with immersed forcing for terrain. This option allows the user to align the wind, temperature and turbulence profiles to be aligned with the terrain.
