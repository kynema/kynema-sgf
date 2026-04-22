.. _capabilities:


Capabilities
============

This section documents a non-exhaustive list of current Kynema-SGF
capabilities.

.. tip::

   If your project relies on a capability that is not yet present in
   Kynema-SGF, please create an issue on the code project page.

   Please acknowledge as a publication co-author any developer that
   has significantly contributed to implementing or improving specific
   capability that was used for that publication.

.. tip::

   The capabilities are linked to the relevant input file references
   (keyword `inp`) and documentation (keyword `doc`). Searching for
   those keywords in the `test/test_files` directory will give
   concrete examples of the feature usage.


Methods and models
------------------

* Numerical methods

   * Advection: second order, piecewise parabolic, piecewise linear, WENO, Bell-Dawson-Shubin [:ref:`doc <discretization>`, :ref:`inp <inputs_incflo_advection>`]

   * Diffusion: second order, explicit, Crank-Nicolson, and implicit [:ref:`doc <discretization>`, :ref:`inp <inputs_incflo_diffusion>`]

   * Mesh refinement: static refinement for specified regions [:ref:`inp <inputs_static_refinement>`], adaptive mesh refinement [:ref:`inp <inputs_amr>`] (e.g., field based, curvature, q-criterion, vorticity [:ref:`inp <inputs_tagging>`])

   * Mesh mapping for non-uniform cartesian grids [:ref:`doc <mapping>`, :ref:`inp <inputs_geometry>`]

* Equations systems

   * Incompressible and low Mach formulations of Navier-Stokes :ref:`[doc] <governing_equations>`

   * Temperature

   * Level set

   * Subgrid scale kinetic energy :ref:`[doc] <turbulence>`

   * Specific dissipation rate :ref:`[doc] <turbulence>`

   * Passive scalar

   * Source terms for these PDEs [:ref:`doc <source_terms>`, :ref:`inp <inputs_momentum_sources>`]

* Turbulence modeling

   * Large Eddy Simulation: constant Smagorinsky,  AMD, one equation :math:`k_{sgs}`, Kosovic [:ref:`doc <turbulence>`, :ref:`inp <inputs_turbulence>`]

   * Wall models: log-law, constant stress, Schumann [:ref:`doc <wall_models>`, :ref:`inp <inputs_abl>`], dynamic (wave model) [:ref:`doc <wall_models>`, :ref:`inp <inputs_boundary_conditions>` 

   * Reynolds-Average Navier-Stokes: :math:`k`-:math:`\omega` SST (and IDDES variant) and One-equation TKE model of Axell and Liungman [:ref:`doc <turbulence>`, :ref:`inp <inputs_turbulence>`]

* Transport models

   * Constant transport coefficients [:ref:`inp <inputs_transport>`]

   * Two phase transport (separate coefficients for each material) [:ref:`inp <inputs_transport>`]

Flow physics
------------

* Wind energy physics

   * Atmospheric boundary layer (ABL): various stability states (stable, unstable, neutral), precursor simulations with inflow boundary planes for wind farm simulations, anelastic formulation, mesoscale forcing, geostrophic forcing, Coriolis forcing, Monin-Obukhov similarity theory, gravity forcing, gravity wave damping [:ref:`inp <inputs_abl>`]

   * Actuator turbine representations: Joukowsky disks, uniform disks, actuator line [:ref:`doc <turbine_models>`, :ref:`inp <inputs_actuator>`]

   * Coupling with OpenFAST

   * Coupling with Nalu-Wind for blade resolved simulations

* Multiphase flows [:ref:`doc <multiphase>`]

   * Prescribed flow cases for verification of volume-of-fluid transport: Zalesak disk, vortex patch

   * Prescribed flow cases for verification of momentum equation coupled to volume-of-fluid transport: Zalesak disk scalar vel, vortex patch scalar vel

   * Validation and demonstration cases: sloshing tank, dam break, breaking waves, falling or inertial droplet

   * Methods to initialize volume-of-fluid field from an initial levelset field

   * Monitors conservation of mass and momentum

* Ocean wave forcing (for multiphase flows) [:ref:`inp <inputs_ocean_waves>`]

   * Wave types: linear (monochromatic), Stokes (second to fifth order), irregular (input by modes files from HOS-Ocean or HOS-NWT)

   * Relaxation zones force wave profile to generate waves at lower x boundary or force toward quiescent flat interface at upper x boundary. Wave profile can also be enforced (instead of numerical beach) at upper x boundary for periodic simulations.

* Boundary conditions

   * Periodic, outflow, inflow, walls, user-defined inflows [:ref:`inp <inputs_boundary_conditions>`]

   * Wall models (e.g., wall functions, stress) [:ref:`doc <wall_models>`, :ref:`inp <inputs_abl>`]

   * Inflow planes from precursor simulations [:ref:`doc <amrwind-abl-bndry-io>`, :ref:`inp <inputs_abl>`]

   * Mesoscale forcing [:ref:`doc <mesoscale_forcing>`, :ref:`inp <inputs_meso_forcing>`]

   * Synthetic turbulence [:ref:`inp <inputs_synthetic_turbulence>`]

   * Arbitrarily spatially and time varying boundary conditions using Python tools [:ref:`inp <inputs_native_boundary_plane>`]

   * Sponge layer driven ABL simulations 

* Geometry

   * Immersed boundary forcing method with stair-case type terrain and wall function [:ref:`doc <terrain>`]

   * Coupling with Nalu-Wind for body-conforming meshes with overset methodology

* Miscellaneous cases

  * Verification and validation cases: method of manufactured solutions, convecting Taylor-Vortex, Rayleigh-Taylor, passive scalar, Burggraf flow, channel flow, Ekman spiral, vortex dipole, vortex ring

* Postprocessing

   * Visualization with VisIt, Paraview, yt

   * Sampling of fields with planes, point probes, lines, volumes, lidar, and radar [:ref:`doc <post_processing>`, :ref:`inp <inputs_sampling>`]

   * Sampling of fields at probes that follow free surface of liquid-gas flows [:ref:`inp <inputs_sampling_freesurface_sampler>`]

   * Scalar outputs such as kinetic energy, enstrophy, total wave energy, and norms [:ref:`doc <post_processing>`, :ref:`inp <inputs_sampling>`]

   * Turbulence averaging quantities such as Reynolds stresses [:ref:`inp <inputs_averaging>`]

   * Field plane averaging and second and third order moments

   * Derived fields and field operators such as vorticity, q-criterion, strain-rates, gradients, divergence, Laplacian [:ref:`inp <inputs_io_derived>`]

   * In-situ post-processing with Ascent

   * Compression for field files [:ref:`doc <compression>`, :ref:`inp <inputs_io>`]

High performance computing
--------------------------

* Highly parallelized and performance portable

   * Shared memory parallelism with OpenMP threading

   * Distributed memory parallelism with MPI

   * Supports all major compilers (e.g., GCC, Intel, LLVM)

   * Runs on all major GPU vendors (NVIDIA, AMD, Intel)

   * Supported build systems: cmake, spack

* Supported linear solvers

   * native AMReX solvers such as MLMG [:ref:`inp <inputs_mlmg>`]

   * hypre

