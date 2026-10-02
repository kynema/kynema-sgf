.. _faq:

Frequently Asked Questions (FAQ) and Common Errors
==================================================

1. "Invalid Particle" error stops the simulation

   This error can appear when running actuator-based turbine simulations coupled to OpenFAST.
   Kynema-SGF uses particles to keep track of actuator points, and the positions of these points
   are determined by the turbine model. If the velocity field in Kynema-SGF becomes
   huge or unrealistic (e.g., not a number), then the turbine model can return bad particle positions
   that do not fit into the Kynema-SGF domain. If this error stops your simulation,
   check the velocity field leading up to the error to determine if something went wrong there.
   If the velocity field is faulty, reevaluate your problem setup.

2. The simulation will not start and the code outputs many lines from the input file

   This happens when a required input file argument is missing. Because AMReX outputs this
   error from every process, and the error message lists all of the lines found in the input
   file, it can be a lot of text to sort through for bigger cases.
   Search for the string "ParmParse" within the output text to locate the expected input
   argument that is missing from the input file so that it can be added.

3. There are many different source terms related to buoyancy or gravity. What do they do and when should they be used?

   The basic buoyancy term often included in ``ICNS.source_terms`` is ``BoussinesqBuoyancy``.
   This source term calculates buoyancy based on the difference between the local temperature
   and a reference temperature, which can be a field or a uniform constant, based on the physics
   setup of the simulation. The resulting buoyancy term is proportional to the thermal expansion
   coefficient (:math:`\beta`).

   When an ABL simulation that has ``BoussinesqBuoyancy`` as a source term has a ``pressure_outflow`` boundary
   condition, the source term ``ABLMeanBoussinesq`` should also be added. This source term helps
   offset the pressure field using a mean temperature profile to be compatible with the pressure boundary
   conditions. If this source term is omitted, nonphysical behavior can appear as the flow nears the pressure outflow
   boundary, especially where thermal stratification occurs.

   When the ``MultiPhase`` physics module is active, the ``GravityForcing`` source term should be
   used. This source term calculates the force of gravity based on the local density. Because temperature
   variations do not cause significant density variations in water, this source term is essential
   for capturing the effect of gravity in a multiphase flow. Adding the input line ``ICNS.use_perturb_pressure = true``
   formulates ``GravityForcing`` as a perturbation from a reference density field. This makes the
   source term compatible with ``pressure_outflow`` conditions as well as other perturbational buoyancy terms.
   For example, a multiphase ABL case would need ``BoussinesqBuoyancy`` activated to include the
   buoyancy effects of temperature in the air, would potentially also need ``ABLMeanBoussinesq`` activated
   if a pressure outflow is being used, and would also need ``GravityForcing`` to capture gravity effect on the
   water with perturbational pressure turned on.