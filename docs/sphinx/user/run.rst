.. _run:

Running Kynema-SGF
=================

To run :program:`kynema_sgf`, the user must provide a text file containing inputs
describing the problem and any additional command-line arguments that override
the parameters in the input file for that particular invocation of the
executable.

.. code-block:: console

   # Parse input parameters from `inputs.abl` but change max_step to 20
   $ ./kynema_sgf inputs.abl time.max_step=20


See :ref:`input file description <inputs>` for a list of Kynema-SGF input flags and proper syntax.

Restarting Kynema-SGF
--------------------

To restart a case, an initial simulation must be run with Kynema-SGF checkpoints. These checkpoints 
are a slice of the solution, and solution settings, that will provide initial input for a restart. Check the
:ref:`time inputs <inputs_time>` and :ref:`i/o inputs <inputs_io>` for input flags to set up checkpoints.

#. Change your ``time.start_time`` to the restart time (optional)
#. Set ``io.restart_file`` to your Kynema-SGF checkpoint directory
#. Re-submit case as above

Restarting Kynema-SGF with OpenFAST Turbines
-------------------------------------------

.. note::
   Kynema-SGF supports OpenFAST version 3.5 and 4.1 (and the patched versions of these). Updates to OpenFAST input files should follow the `OpenFAST documentation <https://openfast.readthedocs.io/en/dev/source/user/api_change.html>`_ on the matter.

.. note::
   Currently, Kynema-SGF automatically creates an OpenFAST .chkp file for every Kynema-SGF checkpoint. These .chkp files 
   are found within each OpenFAST turbine directory.

#. Change your ``time.start_time`` to the restart time (optional)
#. Set ``io.restart_file`` to your Kynema-SGF checkpoint directory
#. Set ``Actuator.T1.openfast_sim_mode = restart`` for each turbine (i.e. T1 in this example)
#. Set ``Actuator.T1.openfast_restart_file`` as the OpenFAST checkpoint file. Note that this must be a relative path from your Kynema-SGF case root, and the ".chkp" must be removed from the filename
#. Set ``Actuator.T1.openfast_start_time`` to the restart time and double check that ``Actuator.T1.openfast_stop_time`` is ok
#. Re-submit case as above
