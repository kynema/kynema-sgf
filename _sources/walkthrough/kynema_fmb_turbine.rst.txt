.. _kynema-fmb-turbine:

Kynema-FMB turbine walkthrough
==============================

This walkthrough shows how to run a turbine simulation in Kynema-SGF with the structural and aerodynamic turbine model provided by Kynema-FMB. The example is based on the regression test input file ``act_kynema_fmb_alm.inp`` and uses the actuator line representation ``TurbineKynemaFMBLine``. To run this example, Kynema-SGF must be compiled with Kynema-FMB support enabled.

Here is the main input file:

.. literalinclude:: ./kynema_fmb_turbines_inp.txt
   :linenos:

* For simplicity, the flow setup here is a simple freestream problem with ``FreeStream Actuator`` enabled in ``incflo.physics``. There is no precursor inflow, wall-model data, or body-force forcing to prepare ahead of time.
* ``Actuator.type`` is set to ``TurbineKynemaFMBLine``, which couples the Kynema-SGF actuator forcing to a Kynema-FMB turbine model.
* The turbine geometry and structural properties are not described directly in the Kynema-SGF input file. Instead, they are read from the WindIO YAML file specified by ``Actuator.TurbineKynemaFMBLine.kynema_fmb_input_file``.
* The actuator discretization still needs to be provided in the Kynema-SGF input file. In this example that includes the blade and tower structural node counts, the number of aerodynamic force points, the Gaussian widths ``epsilon`` and ``epsilon_tower``, and the turbine base locations.
* The turbine solver time step ``Actuator.TurbineKynemaFMBLine.dt`` must divide the Kynema-SGF time step exactly. In this case ``time.fixed_dt = 0.02`` and ``dt = 0.005``, so Kynema-FMB advances four sub-steps per flow solve. Kynema-FMB is highly robust, and it typically can use the same time step size as the flow solver.

The turbine-specific lines in the example are:

* ``Actuator.TurbineKynemaFMBLine.kynema_fmb_input_file = NREL-5MW-aero.yaml`` points to the WindIO description consumed by Kynema-FMB. This file can be found in the regression test directory under ``tests/test_files/actuator_kynema_fmb_alm``; it is not listed here for the sake of brevity.
* ``num_struct_nodes_blade`` and ``num_struct_nodes_tower`` set the number of structural nodes used by the beam model.
* ``num_points_blade`` and ``num_points_tower`` set the aerodynamic sampling points seen by Kynema-SGF. These must match the number of aerodynamic sections available in the Kynema-FMB input data.
* ``rot_speed_rpm`` is for turbine initialization in Kynema-FMB side
* ``density`` is the fluid density, which is needed for the aerodynamic force calculations in Kynema-FMB. This should match the density used in the Kynema-SGF simulation.
* Optional Kynema-FMB solver controls are read from the ``KynemaFMB`` namespace, including ``damping_factor``, ``max_nonlinear_iterations``, ``abs_err_tol``, and ``rel_err_tol``.
* ``Actuator.labels`` together with the per-turbine ``base_position`` entries determine how many turbines are instantiated and where they are placed.

Unlike the OpenFAST walkthrough, this example does not require a separate OpenFAST case directory for each turbine. The run directory instead needs the Kynema-SGF input file and the WindIO YAML file referenced by ``kynema_fmb_input_file``. For the test case shown here, the essential files are:

.. code-block:: console

    act_kynema_fmb_alm.inp  NREL-5MW-aero.yaml

When the run starts, Kynema-SGF builds the flow problem, instantiates one Kynema-FMB turbine for each actuator label, and then advances the two solvers together. Kynema-FMB writes its own turbine output files using names of the form ``kynema_fmb_<label>`` and also writes per-turbine checkpoints whenever Kynema-SGF writes a checkpoint.

Submit the coupled simulation in the usual way:

.. code-block:: console

    kynema_sgf act_kynema_fmb_alm.inp

which should follow after ``srun``, ``mpiexec``, or a similar launcher when running in parallel.

Restarting a Kynema-FMB coupled run
-----------------------------------

Restarting the coupled case requires both the Kynema-SGF checkpoint and the Kynema-FMB turbine checkpoints. The regression test restart file shows the minimal changes:

.. literalinclude:: ./kynema_fmb_turbines_restart_inp.txt
   :linenos:

The restart input works as follows:

* ``FILE = ../act_kynema_fmb_alm/act_kynema_fmb_alm.inp`` reuses the original setup and overrides only the restart-specific entries.
* ``io.restart_file = ../act_kynema_fmb_alm/chk00003`` tells Kynema-SGF to restart the flow solution from checkpoint step 3.
* ``Actuator.TurbineKynemaFMBLine.kynema_fmb_restart_directory`` points to the directory that contains the Kynema-FMB turbine checkpoints.
* ``Actuator.TurbineKynemaFMBLine.kynema_fmb_restart_step = 3`` selects which turbine checkpoint to load.

Internally, Kynema-SGF constructs a Kynema-FMB checkpoint filename for each turbine label using the pattern ``kynema_fmb_00003_<label>.chk``. For the two-turbine example above, restarting from step 3 requires files such as:

.. code-block:: console

    chk00003
    kynema_fmb_00003_T0.chk
    kynema_fmb_00003_T1.chk

all available under the paths referenced by the restart input.

The restart run is then launched the same way as any other case:

.. code-block:: console

    kynema_sgf act_kynema_fmb_alm_restart.inp

If either the AMR checkpoint or one of the turbine checkpoint files is missing, initialization will fail before the simulation begins. 

For additional parameter details beyond this example, see the actuator input reference in the user documentation together with the Kynema-FMB input names used in this walkthrough.