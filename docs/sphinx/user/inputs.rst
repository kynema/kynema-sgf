.. _inputs:

AMR-Wind inputs file
=====================

To run :program:`amr_wind`, the user must provide a text file containing inputs
describing the problem and any additional command-line arguments that override
the parameters in the input file for that particular invocation of the
executable.

.. code-block:: console

   # Parse input parameters from `inputs.abl` but change max_step to 20
   $ ./amr_wind inputs.abl time.max_step=20


The input file is a simple text file containing ``key = value``
entries for the input parameters. The text file can include comments,
any text beginning with ``#`` till the end of line (EOF) is
interpreted as comments and ignored by the parser. Input file
processing is handled by `AMReX ParmParse library
<https://amrex-codes.github.io/amrex/docs_html/Basics.html#parmparse>`_. This
section documents the common input file parameters and their default
values (if available). The focus is on AMR-Wind specific input options
and some common AMReX options (these are fully documented in `AMReX
<https://amrex-codes.github.io/amrex/docs_html/RuntimeParameters.html>`_).

Note that AMReX input files may also include the ``FILE = /path/to/other_file.inp``
directive, which allows users to combine multiple input files. Adding this
directive to an input file (or on the command line) is equivalent to
copying and pasting the contents of ``other_file.inp`` onto that line of the
input file. Unlike standard input file keywords, the directive is applied each
time it is included, rather than only the last having an effect.

In :program:`amr_wind`, the input file
is broken into *sections* indicated by a namespace prefix. For
example, all inputs related to the problem domain are prefixed with
``geometry.`` and so on. A sample input file is shown below

.. literalinclude:: ./amr_wind_inputs.txt
   :linenos:

.. _input-file-ref:

Input file reference
---------------------

The AMR-Wind input file is organized in the following sections

======================= ============================================================
Section                 Description
======================= ============================================================
``geometry``            Computational domain information
``amr``                 Mesh refinement controls
``tagging``             Static and dynamic refinement options
``time``                Simulation time controls
``io``                  Input/Output controls
``incflo``              CFD algorithm and physics controls
``transport``           Transport equation controls
``turbulence``          Turbulence model controls
``ABL``                 Atmospheric boundary layer (ABL) controls
``ABLMesoForcing``      Mesoscale ABL forcing controls
``SyntheticTurbulence`` Inject turbulence using body forces
``Momentum sources``    Activate Momentum source terms and their parameters
``Boundary conditions`` Boundary condition types and gradients
``MLMG options``        Multi-Level Multi-Grid Linear solver options
``Sampling``            Data probes to sample field data during simulations
``Averaging``           Time averaging and correlations
======================= ============================================================

This section documents the parameters available within each section. Please
note that the documentation provided here is for the latest major release of
AMR-Wind. While input file specifications rarely change, major releases of
AMR-Wind (e.g ``2.x`` to ``3.x``) might have breaking changes and the
documentation provided here might not work with older releases.

.. note::

   - Boolean flags (``true/false``) can also be indicated using integers in the text file
     and uses the convention ``0 = False`` and ``1 = True``.

   - Quotes around strings are optional.

   - If an input is repeated only the last one is used.


.. toctree::
   :maxdepth: 2

   inputs_geometry.rst
   inputs_amr.rst
   inputs_tagging.rst
   inputs_time.rst
   inputs_io.rst
   inputs_incflo.rst
   inputs_transport.rst
   inputs_turbulence.rst
   inputs_Momentum_Sources.rst
   inputs_Temperature_Sources.rst
   inputs_TKE_Sources.rst
   inputs_ABL.rst
   inputs_ABL_meso_forcing.rst
   inputs_SyntheticTurbulence.rst
   inputs_Boundary_conditions.rst
   inputs_MLMG.rst
   inputs_post_processing.rst
   inputs_Sampling.rst
   inputs_Subvolume.rst
   inputs_Averaging.rst
   inputs_KineticEnergy.rst
   inputs_Enstrophy.rst
   inputs_FieldNorms.rst
   inputs_Actuator.rst
   inputs_multiphase.rst
   inputs_ocean_waves.rst
