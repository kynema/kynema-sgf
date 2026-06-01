.. _inputs_field_boundaries:
   
Section: Field Boundaries
~~~~~~~~~~~~~~~~~~~~~~~~~
   
Note that field boundaries must be listed to be activated, as described in :input_param:`incflo.field_boundaries`.
There are exceptions to this, however, due to backwards compatibility (BoundaryPlane and ModulatedPowerLaw field
boundaries can be activated through legacy ABL input arguments) or necessity (OceanWavesBoundary is automatically
activated when OceanWaves physics is active).

**BoundaryPlane**

This type of field boundary allows boundary conditions to be specified according to field data contained in planes
that are read in by the solver. These planes are typically created by a precursor simulation, which is usually periodic,
in which case the boundary plane utility does not function as a boundary condition but as a plane data output utility. 
Alternatively, the planes can instead be written by tools external to the solver, allowing data from other sources to
be used to specify boundary conditions. This field boundary does not assume any particular fields for its application
but instead allows the user to specify variable names.

.. input_param:: BoundaryPlane.io_mode 

   **type:** Int, mandatory
   
   Mode for input/output of BoundaryPlane field boundary. A value of 0 indicates that boundary planes should
   be written to file(s) during the simulation, and a value of 1 indicates that boundary planes should be read
   from file(s) during the simulation.

.. input_param:: BoundaryPlane.file

   **type:** String, mandatory

   File name for input/output of BoundaryPlane field boundary. For the native
   output format, this file represents a directory containing plane data, whereas for the netcdf output format,
   this file is a single netcdf file containing all plane data.

.. input_param:: BoundaryPlane.var_names
   
   **type:** List of strings, mandatory

   List of variable names to read/write for the BoundaryPlane field boundary. The corresponding variables must
   be present in the simulation. In read mode, the BoundaryPlane BC is limited to these variables;
   variables not included in this list will use other BCs according to the input file.

.. input_param:: BoundaryPlane.planes

   **type:** List of strings, mandatory

   List of plane names, identifying domain boundaries, to read/write for the BoundaryPlane field boundary. The plane
   names must be in the format "xlo", "xhi", "ylo", "yhi", "zlo", or "zhi". 

.. input_param:: BoundaryPlane.output_format

   **type:** String, optional, default = "native"

   Input/output format for BoundaryPlane field boundary. The default "native" format outputs plane data in a directory
   containing files for each variable and plane. This format corresponds to the AMReX plotfile format, enabling direct 
   visualization through third-party tools like ParaView. The "netcdf" format outputs plane data in a single netcdf file.

.. input_param:: BoundaryPlane.write_frequency

   **type:** Integer, optional, default = 1

   Frequency (actually a time step interval) for writing BoundaryPlane data to file in write mode. This input parameter
   is only relevant to write mode.

.. input_param:: BoundaryPlane.output_start_time

   **type:** Real, optional, default = 0.0

   Time at which to start writing BoundaryPlane data to file in write mode.

.. input_param:: BoundaryPlane.output_and_use_initial_plane

   **type:** Boolean, optional, default = false

   This flag controls whether to output the initial plane at the start of the simulation and use it as the
   boundary plane for the duration of the simulation. If true, this flag implies that the boundary plane is static.

.. input_param:: BoundaryPlane.is_static

   **type:** Boolean, optional, default = false

   This flag controls whether the boundary plane is static or dynamic. If true, the same plane will be used
   for the duration of the simulation. If false, new planes will be read in as time progresses, updating the boundary
   conditions according to the available data. Note that this input parameter is only relevant to read mode.

.. input_param:: BoundaryPlane.fluid_phase

   **type:** String, optional, default = "both"

   This input parameter controls which fluid phase(s) the BoundaryPlane field boundary is applied to when in read mode.
   Valid options are "liquid" (or "water" or "1"), "gas" (or "air" or "2"), and "both" (or "agnostic" or "0").
   This capability works by checking the vof value within the boundary and only using parts of the boundary plane where
   the specified phase is present. As a result, the vof variable must be present for this capability, which 
   means the MultiPhase physics class must be active.

**ModulatedPowerLaw**

This type of field boundary applies a modulated power law profile for velocity
and a linear profile for temperature. The velocity profile can change direction
according to a constant rate, and the temperature profile can incorporate randomness
via a Gaussian distribution. This field boundary is typically used for ABL inflow
conditions, but it can be used independently of the ABL physics class. However,
it does require that velocity and temperature fields are present in the simulation
and automatically applies to all inflow boundary conditions of those fields,
and it does rely on the input arguments :input_param:`ABL.temperature_heights`
and :input_param:`ABL.temperature_values` to specify the reference linear temperature profile.

.. input_param:: MPL.wind_speed

   **type:** Real, optional, default = 8.0

   Reference wind speed for ModulatedPowerLaw field boundary.

.. input_param:: MPL.wind_direction
   
   **type:** Real, optional, default = 270.0

   Wind direction in degrees for ModulatedPowerLaw field boundary. The angle is measured clockwise from the north, so a value of 270 corresponds to a westerly wind.

.. input_param:: MPL.zref

   **type:** Real, optional, default = 90.0

   Reference height for ModulatedPowerLaw field boundary.

.. input_param:: MPL.shear_exp

   **type:** Real, optional, default = 0.1

   Shear exponent for ModulatedPowerLaw field boundary.

.. input_param:: MPL.umax_factor
   
   **type:** Real, optional, default = 1.2

   Maximum velocity factor for ModulatedPowerLaw field boundary.

.. input_param:: MPL.bulk_velocity_factor

   **type:** Real, optional, default = 15.0

   Bulk velocity variable for ModulatedPowerLaw field boundary.

.. input_param:: MPL.shearlayer_height

   **type:** Real, optional, default = 600.0

   Shear layer height for ModulatedPowerLaw field boundary.

.. input_param:: MPL.shearlayer_smear_thickness

   **type:** Real, optional, default = 30.0

   Shear layer smear thickness for ModulatedPowerLaw field boundary.

.. input_param:: MPL.degrees_per_second

   **type:** Real, optional, default = 0.02

   Rate at which to change the wind direction in degrees per second for ModulatedPowerLaw field boundary.

.. input_param:: MPL.start_time

   **type:** Real, optional, default = 0.0

   Time at which to start applying the rate of change in wind direction for ModulatedPowerLaw field boundary.

.. input_param:: MPL.end_time

   **type:** Real, optional, default = very large number

   Time at which to end applying the rate of change in wind direction for ModulatedPowerLaw field boundary.

.. input_param:: MPL.delta_t

   **type:** Real, optional, default = 0.8

   Time scale to use when applying randomness to the temperature inflow of ModulatedPowerLaw field boundary.

.. input_param:: MPL.theta_cutoff_height

   **type:** Real, optional, default = 250.0

   Height below which to apply randomness to the temperature inflow of ModulatedPowerLaw field boundary.

.. input_param:: MPL.theta_gauss_mean

   **type:** Real, optional, default = 0.0

   Mean of the Gaussian distribution to use when applying randomness to the temperature inflow of ModulatedPowerLaw field boundary.

.. input_param:: MPL.theta_gauss_var

   **type:** Real, optional, default = 1.0

   Variance of the Gaussian distribution to use when applying randomness to the temperature inflow of ModulatedPowerLaw field boundary.

**OceanWavesBoundary**

This type of field boundary applies wave profiles for velocity,
volume fraction, and density according to the :ref:`OceanWaves <inputs_ocean_waves>` physics class.
This field is automatically activated when OceanWaves physics is active,
and it relies on fields generated by that physics class, so no specific
input parameters are attributed to this field boundary. In multiphase
simulations, i.e., when the :ref:`MultiPhase <inputs_multiphase>` physics class is active, the
OceanWavesBoundary field boundary applies to all inflow or wave generation boundaries.
At these boundaries, the vof field is populated according to the target wave profile,
and the density field is populated according to the vof field and the specified densities of the two phases.
The velocity field is populated according to the target wave profile as well, but only within the liquid phase.
Whatever other boundary conditions in the input file remain active in the gas phase.
In single phase simulations, i.e., when the MultiPhase physics class is not active and
the waves are intended to be represented as moving terrain, there is no vof field and the density field does
not need to be modified. Therefore, only the velocity within the waves is populated in the boundary condition.
Finally, to avoid conflicts between the OceanWavesBoundary field boundary and the BoundaryPlane field boundary,
the presence of the latter automatically deactivates the former.