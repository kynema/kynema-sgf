.. _inputs_field_boundaries:
   
Section: Field Boundaries
~~~~~~~~~~~~~~~~~~~~~~~~
   
Note that field boundaries must be listed to be activated, as desccribed in :input_param:`incflo.field_boundaries`.
There are exceptions to this, however, due to backwards compatibility (BoundaryPlane and ModulatedPowerLaw field
boundaries can be activated through legacy ABL input arguments) or necessity (OceanWavesBoundary is automatically
activated when OceanWaves physics is active).

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

.. input_param:: BoundaryPlane.output_and_use_initial_plane

   **type:** Boolean, optional, default = false

   This flag controls whether to output the initial plane at the start of the simulation and use it as the
   boundary plane for the duration of the simulation. If true, this flag implies that the boundary plane is static.

.. input_param:: BoundaryPlane.is_static

   **type:** Boolean, optional, default = false

   This flag controls whether the boundary plane is static or dynamic. If true, the same plane will be used
   for the duration of the simulation. If false, new planes will be read in as time progresses, updating the boundary
   conditions according to the available data. Note that this input parameter is only relevant to read mode.

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