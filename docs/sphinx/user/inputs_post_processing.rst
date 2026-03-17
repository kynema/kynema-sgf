.. _inputs_post_processing:

Section: Post-processing
~~~~~~~~~~~~~~~~~~~~~~~~

This section controls post-processing routines supported within
AMR-wind, which include Sampling, Subvolume, Reynolds Averaging (ReAveraging),
ReynoldsStress, TimeAveraging, Enstrophy, FieldNorms, KineticEnergy, and WaveEnergy.

Note that while the input parameters use the keyword ``postproc``, the
actual keyword is determined by the labels provided to
:input_param:`incflo.post_processing`. So for example, if
``incflo.post_processing = my_postproc``, then the options must be prefixed with
``my_postproc.``.

For post-processing If plotfile output is active,
   a plotfile will be written at the end of a simulation (when
   :input_param:`time.stop_time` or :input_param:`time.max_step` is reached), regardless
   of the output timing parameters.

.. input_param:: postproc.type

   **type:** String, optional, default = Sampling

   Specify the type of post-processing routine to apply to this label.

Averaging types (ReAveraging, ReynoldsStress, and TimeAveraging) provide
data via fields that are output to plotfiles; therefore, these types do
not have specific outputs to files in the post_processing directory. However,
for the rest of these routines, which do output data to post-processing files,
there are general input arguments to designate when to write these files,
and these arguments are described below. Note that for these post-processing types,
output will automatically occur at the end of a simulation (when
:input_param:`time.stop_time` or :input_param:`time.max_step` is reached), regardless
of the output timing parameters.

.. input_param:: postproc.output_interval

   **type:** Integer, optional, default = 10

   Specify the output interval (in time steps) when post-processing is performed
   and output to disk. This quantity can instead be specified as ``output_frequency``,
   which was the legacy input argument. Because ``output_interval`` is a more
   accurate name, this is the preferred input argument. 

.. input_param:: postproc.output_time_interval

   **type:** Real, optional

   Specify the output interval (in seconds of simulation time) when post-processing is
   performed and output to disk. The output interval must be specified either in
   simulation time (using this argument) or in time steps.

.. input_param:: postproc.output_delay

   **type:** Integer, optional, default = 0

   Specify the output delay (in time steps) when post-processing and output will begin
   during a simulation. E.g., a delay of 100 will wait until the hundredth timestep to
   check if, according to the output interval, post-processing should be output to disk.

.. input_param:: postproc.output_time_delay

   **type:** Integer, optional, default = 0.

   Specify the output delay (in seconds of simulation time) when post-processing and output will begin
   during a simulation. This parameter is only active in conjunction with an output time interval.

.. input_param:: postproc.enforce_output_time_dt

   **type:** Boolean, optional, default = false

   In the case of a variable dt simulation, the simulation time will not likely correspond
   exactly to the post-processing time interval. However, by setting this parameter to true, the time step size (dt)
   will be shortened when necessary to enforce the simulation time to match the output time interval. 

.. input_param:: postproc.enforce_output_time_dt_reltol

   **type:** Real number, optional, default = 1e-3

   When :input_param:`postproc.enforce_output_time_dt` is true, a tolerance is needed to determine when
   it is necessary to shrink the time step size. This tolerance is relative to the output time interval.
   In most cases, this parameter need not be modified, but it can be changed by the user.

.. input_param:: postproc.output_from_restart

   **type:** Boolean, optional, default = false

   By default, outputs will take place when the total time or time step matches up with
   the specified time or time step interval. However, when restarting a simulation, 
   activating this option allows the relative time to be used, considering the time or time step
   when the current simulation began. This option is not compatible with specifying
   :input_param:`postproc.output_start` or :input_param:`postproc.output_start_time` directly.

.. input_param:: postproc.output_start

   **type:** Integer, optional

   When :input_param:`postproc.output_interval` is active, outputs will take place when the difference
   between the current time step and the initial time step matches up with the specified interval.
   By default, the initial time step used in this calculation is 0. This input argument
   allows the user to override the default behavior by manually specifying
   the initial time step to consider.

.. input_param:: postproc.output_start_time

   **type:** Real number, optional

   When :input_param:`postproc.output_time_interval` is active, outputs will take place when the difference
   between the current time and the initial time matches up with the specified time interval.
   By default, the initial time used in this calculation is 0. This input argument
   allows the user to override the default behavior by manually specifying
   the initial time to consider.

.. input_param:: postproc.output_after_final_step

   **type:** Boolean, optional, default = false

   Similar to checkpoint and plot files, the code can write post-processing
   outputs at the conclusion of a simulation, i.e., after its final step. 
   However, unlike checkpoint and plot files, this behavior is off by default.
   To turn on the final output, the user must specify this input as true.