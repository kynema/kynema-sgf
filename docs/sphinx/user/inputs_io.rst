.. _inputs_io:

Section: io
~~~~~~~~~~~~~~~~~

This section deals with parameters that affect input/output to the simulation,
solely the checkpoint and plot files. These inputs do not affect when these files
are output, but they address file naming and high-level parameters. The "time" section
controls when these files are output.

| Primary location in code: ``kynema-sgf/src/utilities/IOManager.cpp``.

.. input_param:: io.check_file

   **type:** String, optional, default = "chk"

   If :input_param:`time.checkpoint_interval` or :input_param:`time.checkpoint_time_interval` is greater than zero this is the name of the checkpoint
   file appended with the current timestep

.. input_param:: io.plot_file

   **type:** String, optional, default = "plt"

   If :input_param:`time.plot_interval` or :input_param:`time.plot_time_interval` is greater than zero this is the name of the plot
   file appended with the current timestep.

.. input_param:: io.restart_file

   **type:** String, optional, default = ""

   If a string is present Kynema-SGF will restart using the specified file in the string. This is the only argument addressing "input" of data to the simulation instead of "output".

.. input_param:: io.post_processing_directory

   **type:** String, optional, default = "post_processing"

   Name of the directory that will contain post processing output (e.g., sampling).

.. input_param:: io.output_default_variables

   **type:** Boolean, optional, default = true

   Based on what fields are active in a simulation, Kynema-SGF generates a list of variables to output to plotfiles (e.g., velocity, density, and p). If these defaults are not desired, this input argument can be set to false.

.. input_param:: io.allow_missing_restart_fields

   **type:** Boolean, optional, default = true

   When initializing a simulation, Kynema-SGF determines which fields are necessary based on the physics and other details in the input file. If a simulation begins with a restart file, it is possible that the restart file has fewer fields than what the new simulation needs, depending on the input arguments. This argument allows the simulation to continue despite the mismatch. If set to "false", the simulation will abort when necessary fields are missing in the restart file.

.. input_param:: io.outputs

   **type:** List of strings, optional, default = ""

   Add variable names to this input argument to add them to the plotfile output. These must be variables that exist in the simulation and consist of real numbers (not integers).

.. input_param:: io.int_outputs

   **type:** List of strings, optional, default = ""

   Add variable names to this input argument to add them to the plotfile output. These must be variables that exist in the simulation and consist of integers (not real numbers).

.. _inputs_io_derived:

.. input_param:: io.derived_outputs

   **type:** List of strings, optional, default = ""

   Add derived variable names to this input argument
   to add them to the plotfile output. These are derived
   quantities that are functions of real variables that exist
   in the simulation. Currently, the available derived quantity definitions
   that operate on the velocity field are vorticity magnitude
   (``mag_vorticity``), q-criterion (``q_criterion``),
   non-dimensional q-criterion (``q_criterion_nondim``),
   and strain rate magnitude (``mag_strainrate``). Generic
   derived quantity definitions, which operate on fields specified as an argument,
   include the gradient operator (``grad``), the divergence
   operator (``div``), the Laplacian operator (``laplacian``),
   and components (``components``), which isolates the specified
   component of a field.


.. input_param:: io.skip_outputs

   **type:** List of strings, optional, default = ""

   Add variable names to this input argument to omit them from the plotfile output. These refer to variables that are be real numbers, and this is a way to individually omit default output variables.

.. input_param:: io.output_hdf5_plotfile

   **type:** Boolean, optional, default = false

   Flag indicating whether or not to output HDF5 plot files.

.. input_param:: io.hdf5_compression

   **type:** String, optional, default = "None@0"

   String for requesting a particular ZFP compression in HDF5 plot files.

.. input_param:: io.nfiles

   **type:** Int, optional, default = 256

   Number of plot and checkpoint data files per write. If the system's IO prefers fewer or more files, this number can be modified with this option.
