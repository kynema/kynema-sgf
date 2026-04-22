.. _inputs_time:

Section: time
~~~~~~~~~~~~~~~~~

This section deals with parameters that control the time stepping of the simulation. 
This section also addresses the time-dependent nature of checkpoint files, plot files, and regridding. 

| Primary location in code: ``amr-wind/core/SimTime.cpp``.

.. input_param:: time.stop_time

   **type:** Real number

   The time (in seconds) when the simulation should stop. This parameter is
   related to :input_param:`time.max_step` depending on whether the user
   provides a positive or negative value for this input. If the user provides a
   positive value, then the stop criteria is such that the simulation stops when
   it reaches the ``stop_time``. If the user specifies a negative value, then
   the simulation stops when it has completed the prescribed number of steps as
   determined by :input_param:`time.max_step`.

.. input_param:: time.max_step

   **type:** Integer

   The number of timesteps to run before termination. See also
   :input_param:`time.stop_time`.

.. input_param:: time.fixed_dt

   **type:** Real number

   Fixed timestep size (in seconds) used to advance the simulation. If this
   parameter is negative, then :input_param:`time.cfl` is used to determine the
   adaptive timestep during the simulation.
   
.. input_param:: time.initial_dt

   **type:** Real number

   Initial timestep size (in seconds) used to initialize the simulation. 
   Only activated if :input_param:`time.fixed_dt` is negative 
   (signaling CFL controlled time stepping) and if :input_param:`time.initial_dt` is positive.
   This parameter can be useful for starting CFL-controlled simulations like 
   Rayleigh-Taylor flow that initialize with zero velocity.

.. input_param:: time.max_dt

   **type:** Real number

   Maximum timestep size (in seconds) used to limit CFL-controlled time stepping. 
   Only activated if :input_param:`time.max_dt` is positive and fixed time stepping is disabled. 
   This parameter can be useful if CFL-controlled timestepping is desired but there
   are additional constraints to the timestep size considerations.

.. input_param:: time.min_dt

   **type:** Real number

   Minimum timestep size (in seconds) used to limit CFL-controlled time stepping by stopping
   a simulation when the timestep size has become too small. 
   Only activated if :input_param:`time.min_dt` is positive and fixed time stepping is disabled.
   This parameter is more useful in the development context when running simulations that are not
   reliably stable. If a simulation's velocity values dramatically grow, i.e., they blow up,
   CFL-controlled timestepping would make the timestep size unrealistically small, often allowing
   the simulation to continue to run indefinitely with an nonphysical flow field. When this parameter
   is active, the code will abort when the timestep size becomes smaller than the specified minimum,
   preventing the simulation from continuing to run wastefully.

.. input_param:: time.cfl

   **type:** Real number, optional, default = 0.5

   The maximum CFL allowed during the simulation. Adaptive timestepping algorithm is 
   used to ensure that the maximum CFL condition is not violated while using the largest
   allowable timestep to advance the simulation.

.. input_param:: time.init_shrink

   **type:** Real number, optional, default = 0.1

   ``init_shrink`` indicates a reduction factor applied to the timestep size
   during initialization phased used to perform the initial iterations to
   determine pressure. This factor must be a value greater than zero but less
   than or equal to 1.0.

.. input_param:: time.max_dt_growth

   **type:** Real number, optional, default = 0.1

   ``max_dt_growth`` indicates a growth limit factor applied to the timestep size
   whenever it increases during a variable dt simulation. This factor must be a 
   value greater than zero but less than or equal to 1.0. In most cases a value of
   ``0.1`` is sufficient, but if dt is forced to shrink harshly and often due to 
   :input_param:`time.enforce_plot_time_dt` or :input_param:`time.enforce_checkpoint_time_dt` 
   being active, a larger growth factor may be more practical.

.. input_param:: time.regrid_interval

   **type:** Integer, optional, default = -1

   If :input_param:`amr.max_level` is greater than zero, this parameter
   indicates the frequency (in timesteps) at which the mesh is adaptively
   refined based on various user-specified criteria. If this value is negative,
   the mesh is only refined once during initialization and remains constant for
   the rest of the simulation.

.. input_param:: time.plot_interval

   **type:** Integer, optional, default = -1

   If this value is greater than zero, it indicates the frequency (in timesteps)
   at which outputs (plot files) are written to disk. If plotfile output is active,
   a plotfile will be written at the end of a simulation (when
   :input_param:`time.stop_time` or :input_param:`time.max_step` is reached), regardless
   of the output timing parameters.

.. input_param:: time.plot_time_interval

   **type:** Real number, optional, default = -1.0

   If this value is greater than zero, it indicates the frequency (in seconds)
   at which outputs (plot files) are written to disk. This cannot be specified in conjunction with :input_param:`time.plot_interval`. 
   Because the size of time steps (dt) do not automatically correspond to the size of the plot time interval,
   additional parameters are available if desired: :input_param:`time.plot_time_interval_reltol`, 
   :input_param:`time.enforce_plot_time_dt`, and :input_param:`time.enforce_plot_dt_reltol`. These additional parameters
   are only relevant for simulations with a variable time step size (dt).

.. input_param:: time.plot_delay

   **type:** Integer, optional, default = 0

   If :input_param:`time.plot_interval` is greater than zero, then the plot delay specifies how long (in timesteps)
   to wait before writing a plot file. The implementation waits until this threshold is reached to check if the interval 
   allows for a file to be written. For example, if the plot delay is specified to be "1000", and the plot
   interval is "10", then the first plot file written would be at timestep "1000". If the plot delay is "1001" and the
   plot interval is still "10", then the first plot file file would be at timestep "1010".

.. input_param:: time.plot_time_delay

   **type:** Real number, optional, default = 0.0

   If :input_param:`time.plot_time_interval` is greater than zero, then the plot time delay specifies how long (in seconds)
   to wait before writing a plot file. Similar to :input_param:`time.plot_delay`, the implementation waits until this threshold is 
   reached to check if the time interval allows for a file to be written.

.. input_param:: time.checkpoint_interval

   **type:** Integer

   If this value is greater than zero, it indicates the frequency (in timesteps)
   at which checkpoint (restart) files are written to disk. If checkpoint output is active,
   a checkpoint will be written at the end of a simulation (when
   :input_param:`time.stop_time` or :input_param:`time.max_step` is reached), regardless
   of the output timing parameters.

.. input_param:: time.checkpoint_time_interval

   **type:** Real number, optional, default = -1.0

   If this value is greater than zero, it indicates the frequency (in seconds)
   at which checkpoint (restart) files are written to disk. This cannot be specified in conjunction with :input_param:`time.checkpoint_interval`. 
   Because the size of time steps (dt) do not automatically correspond to the size of the checkpoint time interval,
   additional parameters are available if desired: :input_param:`time.checkpoint_time_interval_reltol`, 
   :input_param:`time.enforce_checkpoint_time_dt`, and :input_param:`time.enforce_checkpoint_dt_reltol`. These additional parameters
   are only relevant for simulations with a variable time step size (dt).

.. input_param:: time.checkpoint_delay

   **type:** Integer, optional, default = 0

   If :input_param:`time.checkpoint_interval` is greater than zero, then the checkpoint delay specifies how long (in timesteps)
   to wait before writing a checkpoint file. The implementation waits until this threshold is reached to check if the interval 
   allows for a file to be written. For example, if the checkpoint delay is specified to be "1000", and the checkpoint
   interval is "10", then the first checkpoint file written would be at timestep "1000". If the checkpoint delay is "1001" and the
   checkpoint interval is still "10", then the first checkpoint file would be at timestep "1010".

.. input_param:: time.checkpoint_time_delay

   **type:** Real number, optional, default = 0.0

   If :input_param:`time.checkpoint_time_interval` is greater than zero, then the checkpoint time delay specifies how long (in seconds)
   to wait before writing a checkpoint file. Similar to :input_param:`time.checkpoint_delay`, the implementation waits until this threshold is 
   reached to check if the time interval allows for a file to be written.
   
.. input_param:: time.regrid_start

  **type:** Integer, optional, default = 0; default = start index upon restart

  This user-specified parameter sets the base timestep onward which the mesh is adaptively
  refined.

.. input_param:: time.plot_start

  **type:** Integer, optional, default = 0; default = start index upon restart

  This user-specified parameter sets the base timestep onward which the output (plot files)
  are written to the disk. This parameter is specifically for offsetting the index following a restart.

.. input_param:: time.checkpoint_start

  **type:** Integer, optional, default = 0; default = start index upon restart

  This user-specified parameter sets the base timestep onward which the checkpoint (restart) 
  files are written to the disk. This parameter is specifically for offsetting the index following a restart.

.. input_param:: time.use_force_cfl

   **type:** Boolean, optional, default = true

   If this flag is true then the forces (including the pressure gradient) are included
   in the CFL calculation.

.. input_param:: time.plot_time_interval_reltol

   **type:** Real number, optional, default = 1e-2

   When :input_param:`time.plot_time_interval` is greater than zero, the implementation 
   compares the current simulation time to the specified time interval in order to output
   the time step that meets or just passes the time interval. Because this involves a comparison
   between real numbers, this comparison uses a tolerance, which is relative to the current timestep,
   which helps avoid skipping time steps which are very close, but just shy, of the time interval.
   In most cases, this parameter need not be modified, but it can be changed by the user.
   
.. input_param:: time.enforce_plot_time_dt

   **type:** Boolean, optional, default = false

   In the case of a variable dt simulation, the simulation time will not likely correspond
   exactly to the plot time interval. However, by setting this parameter to true, the time step size (dt)
   will be shortened when necessary to enforce the simulation time to match the plot time interval. Enabling 
   this feature while also enabling :input_param:`time.enforce_checkpoint_time_dt` will result in a warning message.

.. input_param:: time.enforce_plot_dt_reltol

   **type:** Real number, optional, default = 1e-3

   When :input_param:`time.enforce_plot_time_dt` is true, a tolerance is needed to determine when
   it is necessary to shrink the time step size. This tolerance is relative to the plot time interval.
   In most cases, this parameter need not be modified, but it can be changed by the user.

.. input_param:: time.checkpoint_time_interval_reltol

   **type:** Real number, optional, default = 1e-2

   This parameter is active when :input_param:`time.checkpoint_time_interval` is greater than zero,
   and it exists for the same reason as :input_param:`time.plot_time_interval_reltol`.

.. input_param:: time.enforce_checkpoint_time_dt

   **type:** Boolean, optional, default = false
   
   Similar to :input_param:`time.enforce_plot_time_dt`, setting this parameter to true will enforce 
   the simulation time to match the checkpoint time interval in the case of a variable dt simulation. 
   Enabling this feature while also enabling :input_param:`time.enforce_plot_time_dt` will result in 
   a warning message.

.. input_param:: time.enforce_checkpoint_dt_reltol

   **type:** Real number, optional, default = 1e-3

   When :input_param:`time.enforce_checkpoint_time_dt` is true, a tolerance is needed to determine when
   it is necessary to shrink the time step size. This tolerance is relative to the checkpoint time interval.
   In most cases, this parameter need not be modified, but it can be changed by the user.

.. input_param:: time.profiling_interval

   **type:** Integer, optional, default = -1

   If this value is greater than zero, it indicates the frequency (in
   timesteps) at which profiling information is written to the log
   file. This will only work if the code was compiled with
   ``KYNEMA_SGF_ENABLE_TINY_PROFILE`` set to ``ON``.
