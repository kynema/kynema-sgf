.. _inputs_actuator_kynema_fmb:

Actuator-Based Kynema-FMB Coupling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Kynema-SGF supports actuator-based coupling to Kynema-FMB through two actuator types:

* ``TurbineKynemaFMBLine`` for the actuator-line representation
* ``TurbineKynemaFMBDisk`` for the actuator-disk representation

Both types share the same Kynema-FMB turbine parser and therefore accept the same input arguments. Unless noted otherwise, every input documented below with the prefix ``Actuator.TurbineKynemaFMBLine`` is also accepted with the prefix ``Actuator.TurbineKynemaFMBDisk`` when ``Actuator.type = TurbineKynemaFMBDisk``.

These coupled turbine types require a Kynema-SGF build with Kynema-FMB support enabled. They also require ``Actuator`` in ``incflo.physics`` and typically require ``ActuatorForcing`` in ``ICNS.source_terms``.

Unlike the OpenFAST-coupled turbine models, the rotor diameter and hub height are read from the WindIO file referenced by ``kynema_fmb_input_file``. They are not separate Kynema-SGF user inputs for these actuator types.

Example for ``TurbineKynemaFMBLine``::

   incflo.physics = FreeStream Actuator
   ICNS.source_terms = ActuatorForcing
   Actuator.labels = T0 T1
   Actuator.type = TurbineKynemaFMBLine

   Actuator.TurbineKynemaFMBLine.kynema_fmb_input_file = NREL-5MW-aero.yaml
   Actuator.TurbineKynemaFMBLine.num_struct_nodes_blade = 11
   Actuator.TurbineKynemaFMBLine.num_struct_nodes_tower = 11
   Actuator.TurbineKynemaFMBLine.num_points_blade = 64
   Actuator.TurbineKynemaFMBLine.num_points_tower = 7
   Actuator.TurbineKynemaFMBLine.num_blades = 3
   Actuator.TurbineKynemaFMBLine.epsilon = 5.0 5.0 5.0
   Actuator.TurbineKynemaFMBLine.epsilon_tower = 5.0 5.0 5.0
   Actuator.TurbineKynemaFMBLine.fllc = false
   Actuator.TurbineKynemaFMBLine.rot_speed_rpm = 12.1
   Actuator.TurbineKynemaFMBLine.dt = 0.005
   Actuator.TurbineKynemaFMBLine.density = 1.225
   Actuator.TurbineKynemaFMBLine.nacelle_drag_coeff = 1.0
   Actuator.TurbineKynemaFMBLine.nacelle_area = 8.0
   Actuator.TurbineKynemaFMBLine.output_frequency = 1

   Actuator.T0.base_position = 0.0 -250.0 -90.0
   Actuator.T1.base_position = 0.0 250.0 -90.0

   KynemaFMB.abs_err_tol = 1e-6
   KynemaFMB.rel_err_tol = 1e-4

For a disk-based run, set ``Actuator.type = TurbineKynemaFMBDisk`` and replace the parameter prefix accordingly.

If the individual turbines (e.g. ``T0`` and ``T1``) differ in any other details, the input arguments can
be specified using the turbine labels (e.g. ``Actuator.T0.nacelle_area``, ``Actuator.T1.nacelle_area``)
instead of the turbine type prefix, which assigns the value to all turbines of that type.

Per-turbine placement
"""""""""""""""""""""

.. input_param:: Actuator.T0.base_position

   **type:** List of 3 real numbers, required

   Base position of the turbine in the global coordinate system. Replace ``T0`` with each entry listed in ``Actuator.labels``.

Kynema-FMB turbine setup and coupling
"""""""""""""""""""""""""""""""""""""

.. input_param:: Actuator.TurbineKynemaFMBLine.kynema_fmb_input_file

   **type:** String, required

   Path to the Kynema-FMB input file in WindIO format. Kynema-SGF reads the rotor diameter and hub height from this file and forwards the structural, aerodynamic, and mass-property data to Kynema-FMB.

.. input_param:: Actuator.TurbineKynemaFMBLine.num_struct_nodes_blade

   **type:** Int, required

   Number of structural nodes used by Kynema-FMB for each blade.

.. input_param:: Actuator.TurbineKynemaFMBLine.num_struct_nodes_tower

   **type:** Int, required

   Number of structural nodes used by Kynema-FMB for the tower.

.. input_param:: Actuator.TurbineKynemaFMBLine.num_points_blade

   **type:** Int, required

   Number of aerodynamic sections per blade used by the actuator representation. This must match the number of aerodynamic sections available in the Kynema-FMB input.

.. input_param:: Actuator.TurbineKynemaFMBLine.num_points_tower

   **type:** Int, required

   Number of aerodynamic tower sections used by the actuator representation. Set this to ``0`` to disable tower aerodynamics. If nonzero, it must match the tower discretization available in the Kynema-FMB input.

.. input_param:: Actuator.TurbineKynemaFMBLine.num_blades

   **type:** Int, optional, default = 3

   Number of blades in the actuator representation.

.. input_param:: Actuator.TurbineKynemaFMBLine.density

   **type:** Real, required

   Fluid density used to non-dimensionalize the forces exchanged with Kynema-FMB. This should match the density in the Kynema-SGF simulation.

.. input_param:: Actuator.TurbineKynemaFMBLine.rot_speed_rpm

   **type:** Real, optional, default = 0

   Initial rotor speed in RPM. This parameter is ignored when :input_param:`Actuator.TurbineKynemaFMBLine.rot_speed_radps` is present.

.. input_param:: Actuator.TurbineKynemaFMBLine.rot_speed_radps

   **type:** Real, optional, default = 0

   Initial rotor speed in radians per second. If present, it overrides :input_param:`Actuator.TurbineKynemaFMBLine.rot_speed_rpm`.

.. input_param:: Actuator.TurbineKynemaFMBLine.generator_power_init

   **type:** Real, optional, default = 0

   Initial generator power (in W) passed to Kynema-FMB.

.. input_param:: Actuator.TurbineKynemaFMBLine.hub_wind_vector_init

   **type:** List of 3 real numbers, optional, default = 0 0 0

   Initial wind vector seen by the turbine hub. Kynema-SGF converts this vector to a speed magnitude and passes it to Kynema-FMB for the sake of the controller.

.. input_param:: Actuator.TurbineKynemaFMBLine.yaw_deg

   **type:** Real, optional, default = 0

   Initial nacelle yaw angle in degrees. This parameter is ignored when :input_param:`Actuator.TurbineKynemaFMBLine.yaw_rad` is present.

.. input_param:: Actuator.TurbineKynemaFMBLine.yaw_rad

   **type:** Real, optional, default = 0

   Initial nacelle yaw angle in radians. If present, it overrides :input_param:`Actuator.TurbineKynemaFMBLine.yaw_deg`.

.. input_param:: Actuator.TurbineKynemaFMBLine.generator_efficiency

   **type:** Real, optional, default = 1

   Generator efficiency passed to Kynema-FMB.

.. input_param:: Actuator.TurbineKynemaFMBLine.controller_shared_library_path

   **type:** String, optional, default = empty

   Path to the controller shared library passed to Kynema-FMB. This is typically a ROSCO dynamic library.

.. input_param:: Actuator.TurbineKynemaFMBLine.controller_input_file

   **type:** String, optional, default = empty

   Path to the controller input file. When this value is provided, Kynema-SGF creates the Kynema-FMB controller interface and expects the shared library path to be set as well.

.. input_param:: Actuator.TurbineKynemaFMBLine.dt

   **type:** Real, optional, default = same as Kynema-SGF time step

   Kynema-FMB time step size. It must divide the Kynema-SGF time step so that the coupled solver can take an integer number of Kynema-FMB sub-steps per CFD time step.

Actuator forcing and output controls
""""""""""""""""""""""""""""""""""""

.. input_param:: Actuator.TurbineKynemaFMBLine.epsilon

   **type:** List of 3 real numbers, conditionally required

   Gaussian smearing width for the blade force projection. Kynema-SGF requires at least one of ``epsilon`` or ``epsilon_chord`` for turbine actuator simulations. When FLLC is enabled, both ``epsilon`` and ``epsilon_chord`` are required.

.. input_param:: Actuator.TurbineKynemaFMBLine.epsilon_chord

   **type:** List of 3 real numbers, optional

   Non-dimensional ``epsilon / chord`` values used to compute a blade-force smearing width from the local chord. Kynema-SGF uses the maximum of the direct ``epsilon`` value and the chord-based value when both are provided.

.. input_param:: Actuator.TurbineKynemaFMBLine.epsilon_min

   **type:** List of 3 real numbers, optional

   Minimum value allowed when chord-based epsilon values are used.

.. input_param:: Actuator.TurbineKynemaFMBLine.epsilon_tower

   **type:** List of 3 real numbers, optional

   Gaussian smearing width for tower forces.

.. input_param:: Actuator.TurbineKynemaFMBLine.fllc

   **type:** Bool, optional, default = false

   Enables the filtered lifting line correction for the blade loading.

.. input_param:: Actuator.TurbineKynemaFMBLine.fllc_type

   **type:** String, optional, default = ``variable_chord``

   Selects the FLLC formulation. Supported values are ``constant_chord`` and ``variable_chord``.

.. input_param:: Actuator.TurbineKynemaFMBLine.fllc_relaxation_factor

   **type:** Real, optional, default = 0.1

   Relaxation factor applied to the FLLC velocity update.

.. input_param:: Actuator.TurbineKynemaFMBLine.fllc_start_time

   **type:** Real, optional, default = 0

   Simulation time at which the FLLC correction becomes active.

.. input_param:: Actuator.TurbineKynemaFMBLine.fllc_nonuniform

   **type:** Bool, optional, default = true

   Enables the non-uniform radial point distribution used by the FLLC implementation.

.. input_param:: Actuator.TurbineKynemaFMBLine.fllc_epsilon_dr_ratio

   **type:** Real, optional, default = 1

   Target ratio of epsilon to actuator-point spacing used when building the non-uniform FLLC point distribution.

.. input_param:: Actuator.TurbineKynemaFMBLine.nacelle_drag_coeff

   **type:** Real, optional, default = 0

   Nacelle drag coefficient used to model nacelle drag.

.. input_param:: Actuator.TurbineKynemaFMBLine.nacelle_area

   **type:** Real, optional, default = 0

   Frontal area of the nacelle used together with ``nacelle_drag_coeff``.

.. input_param:: Actuator.TurbineKynemaFMBLine.output_frequency

   **type:** Int, optional, default = 10

   Frequency, in Kynema-SGF time steps, for writing the actuator NetCDF output produced by Kynema-SGF. Kynema-FMB's own turbine output is written separately every Kynema-SGF time step.

Restart and controller options
""""""""""""""""""""""""""""""

For a full coupled restart, these turbine-specific settings are used together with the standard Kynema-SGF flow restart input ``io.restart_file``.

.. input_param:: Actuator.TurbineKynemaFMBLine.kynema_fmb_restart_step

   **type:** Int, optional

   External restart step for Kynema-FMB. When this value is present, Kynema-SGF switches the turbine to restart mode and constructs a per-turbine checkpoint filename of the form ``kynema_fmb_00003_<label>.chk``.

.. input_param:: Actuator.TurbineKynemaFMBLine.kynema_fmb_restart_directory

   **type:** String, optional, default = ``.``

   Directory containing the Kynema-FMB checkpoint files referenced by ``kynema_fmb_restart_step``.

Kynema-FMB solver options
"""""""""""""""""""""""""

.. input_param:: KynemaFMB.abs_err_tol

   **type:** Real, optional, default = 1e-5

   Absolute error tolerance used by the Kynema-FMB nonlinear solver.

.. input_param:: KynemaFMB.rel_err_tol

   **type:** Real, optional, default = 1e-4

   Relative error tolerance used by the Kynema-FMB nonlinear solver.

.. input_param:: KynemaFMB.max_nonlinear_iterations

   **type:** Int, optional, default = 12

   Maximum number of nonlinear iterations used by the Kynema-FMB solver per sub-step.

.. input_param:: KynemaFMB.damping_factor

   **type:** Real, optional, default = 0

   Numerical damping factor forwarded to the Kynema-FMB solver, which is applied in its temporal scheme.
   Counterintuitively, 0 corresponds to full damping and 1 corresponds to no damping.