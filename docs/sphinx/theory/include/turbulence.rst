.. _turbulence:

Turbulence Models
-----------------

RANS models
~~~~~~~~~~~~~

The RANS models are available in two flavors: wall-modeled and wall-resolved. The former model is 
designed for cases with :math:`y+ > 30` while the latter requires :math:`y+ < 5`. The wall-modeled RANS 
model available in Kynema-SGF is based on the work of `Axell and Liungman (EFM 2001 ) <https://link.springer.com/article/10.1023/A:1011560202388>`_.
The code also includes Menter's K-Omega SST model with IDDES support. 

Axell One-Equation RANS Model 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The one-equation model solves the transport equation for turbulent kinetic energy (TKE). The length scale is computed using algebraic equations. 
The transport equation for TKE is given by: 

.. math:: \frac{ \partial (\rho k)}{\partial t} 
   + \nabla \cdot (\rho U k)  =  P_s + P_b - \epsilon + D 

Here :math:`P_s` is the shear production term, :math:`P_b` is the buoyancy production/destruction term, :math:`\epsilon` is the turbulent dissipation 
rate and :math:`D` is the turbulent diffusion term. These terms are computed as follows 

.. math:: P_s= \nu_t S^2

.. math:: P_b= -{\nu_t}^{'} N^2 

.. math:: \epsilon= C_0 \frac{k^3/2}{L}

.. math:: D = \frac{\partial}{\partial x_j} [(\nu+\nu_t)\frac{\partial U_j}{\partial x_i}]

Here :math:`P_s` is the strain rate, :math:`P_b` is the buoyancy frequency and :math:`L` is the length scale computed algebraically.
The strain rate and buoyancy frequency are computed using the same method used in the literature and are not repeated here. The 
length scale is computed as follows: 

.. math:: \frac{1}{L^2} = \frac{1}{Ls^2} + \frac{1}{Lb^2}

The shear length scale is given by :math:`Ls=\kappa z`. An upper limit can be imposed for the shear length scale to avoid excessive values. 
In the current model, it is set to 30 and can be modified to be computed from Geostrophic wind too. The buoyancy length scale is given by 

.. math:: Lb = Cb \frac{\sqrt{k}}{N} 

The implementation methodology is different for stable/neutral and unstable stratification and follows the recommendation in the paper. The
turbulent viscosity is computed as follows: 

..  math:: \nu_t = C_\mu \sqrt{k} L 

..  math:: {\nu_t}^{'} = {C_\mu}^{'} \sqrt{k} L 

Here :math:`C_\mu` and :math:`{C_\mu}^{'}` are non-uniform model constants which depend on :math:`C_0` and turbulent Richardson number 
:math:`Rt`. The calculations of these terms can be found in the reference. The turbulent Prandtl number also depends on the turbulent 
Richardson number and is computed using am empirical expression from the reference. The boundary condition for TKE at the lower boundary 
is given by: 

.. math:: k = k_w ^ {(2/3)}

.. math:: k_w = \frac{{u_*}^{3}}{{C_0}^3} + \frac{\max{(Q,0)}\kappa d_1}{{C_0}^3}

Here :math:`Q` is the sensible heat flux at the surface and :math:`d_1`  is the near-wall distance. For cases with terrain, there is also 
a check for near-wall distance from the surface of the terrain. The wall boundary condition is implemented as a forcing term at the first cell 
above the lower surface and terrain. 


Separation treatments for the Axell model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``KLAxellSeparation`` turbulence model is the Axell one-equation model
above with optional treatments for boundary-layer flow that separates over
terrain. With every treatment disabled it gives the same results as
``KLAxell``, and on flat terrain the treatments that are gated by the sensor
below leave the solution unchanged. The treatments are selected in the
``KLAxellSeparation`` input section and tuned in
``KLAxellSeparation_coeffs`` (see :ref:`inputs_turbulence`).

**Sensor and gate.** The pressure-gradient sensor

.. math::

   s = \frac{\hat{u}_i \, \partial p / \partial x_i}
            {\rho \, (k + c_u |\mathbf{u}|^2) / L}

marks cells where the flow runs against an adverse pressure gradient. The
treatments that depend on it act through a gate :math:`g` that ramps from 0 at
the threshold :math:`s_T` to 1 at :math:`2 s_T`. By default the gate relaxes
toward that ramp value with a time scale of 10 s
(``KLAxellSeparation_coeffs.gate_relaxation_time``). An instantaneous gate lets
the treatments, the pressure field and the sensor feed back on each other
from one time step to the next, which produced grid-scale stripes in the eddy
viscosity. A gate that is 0 stays exactly 0 while the sensor stays below the
threshold. Once the relaxed gate has opened, it decays toward 0 with the time
scale :math:`\tau` after the sensor falls below the threshold, so it stays
positive for a while.

**Treatments.**

- ``realizable_cmu``: the eddy viscosity and the shear and buoyancy
  production are divided by :math:`1 + c_s \, g \, \max(0, \Sigma / C_\mu - 1)`,
  with :math:`\Sigma = L S / \sqrt{k}`, which equals :math:`C_\mu` in an
  equilibrium log layer. This limits the eddy viscosity where the strain rate
  is large compared with the turbulence, as in an adverse pressure gradient.
- ``production_cap``: the shear production above :math:`C_P \varepsilon` is
  removed where the gate is open. It acts as a safety limit at stagnation
  points.
- ``destruction_boost``: the dissipation in the source of the turbulent
  kinetic energy is multiplied by :math:`1 + g (c_d - 1)`.
- ``curvature_correction``: a streamline-curvature correction that is not
  gated by the sensor and does nothing where the streamlines are straight. The
  default ``richardson`` model adds a curvature Richardson number to
  :math:`R_t`; the ``rotation_function`` model is also available.
- ``implicit_dissipation``: a numerical option that treats the dissipation
  implicitly in the diffusion solve. Converged solutions are unchanged; it
  prevents the dissipation from driving the turbulent kinetic energy negative
  with large time steps.

**Recommended settings.** For flow over terrain where separation matters:

.. list-table::
   :header-rows: 1

   * - Input
     - Recommendation
   * - ``KLAxellSeparation.pressure_gradient_sensor`` and
       ``KLAxellSeparation.realizable_cmu``
     - ``true``; the main treatment.
   * - ``KLAxellSeparation_coeffs.gate_relaxation_time``
     - Keep the default of 10 s. A value of 0 gives the instantaneous gate.
   * - ``KLAxellSeparation.curvature_correction``
     - Optional; ``true`` adds the effect of streamline curvature, which also
       changes the flow away from the separation region.
   * - ``KLAxellSeparation.destruction_boost``,
       ``KLAxellSeparation.production_cap``
     - Optional; small effect with the default coefficients.
   * - ``KLAxellSeparation.implicit_dissipation``
     - Optional, for large time steps.
   * - ``KLAxellSeparation.sensor_source``
     - Keep ``pressure``. The ``velocity`` sensor concentrated the limiter in
       a thin band at crest height.

**Smooth-hill test.** A periodic two-dimensional hill of height 100 m and
half-width 400 m in a neutral boundary layer (``TerrainDrag``, 4 m grid,
300 s) separates on the lee slope with ``KLAxell``. The number of cells with
reversed flow was 117 with ``KLAxell``, 135 with the sensor, the limiter and
the relaxed gate, 157 with the curvature correction added, and 165 with all
treatments on. The instantaneous gate gave the same separation as the relaxed
gate but a cell-to-cell roughness of the eddy viscosity about 190 times that
of ``KLAxell``; the relaxed gate reduced it to about 5 times, as one smooth
region under the shear layer that leaves the crest. The same trend held on an
8 m grid (11, 20 and 33 cells). Away from the hill the eddy viscosity changed
by less than about 10 percent for most cells. These runs show the effect of the
treatments; they are not a comparison with measurements.

The shear production passed to the production cap carries the density of the
eddy viscosity, so :math:`C_P` is measured against :math:`\rho \varepsilon`.

LES models for subgrid scales
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Smagorinsky model
^^^^^^^^^^^^^^^^^

Simple eddy viscosity model, the dissipation is calculated using the
resolved strain rate tensor and the grid resolution as

.. math::

   \begin{aligned}
       \tau_{ij} &= -2 \nu_t \widetilde{S}_{ij} \\
       \nu_t &= C_s^2 \Delta^2 (2 \langle S_{ij} S_{ij} \rangle)^{\frac{1}{2}}
   \end{aligned}


AMDNoTherm model
^^^^^^^^^^^^^^^^^
This is the implementation of the base AMD model, useful for flows without a temperature field.

The eddy viscosity is calculated using an anisotropic derivative with a
different filter width in each direction

.. math::

   \begin{aligned}
       \hat{\partial}_i &= \sqrt{C} \delta_i \partial_i \textrm{ for } i=1,2,3 \\
       C &= 1/3, \textrm{ Poincare coefficient for } 2^{nd} \textrm{ order gradient} \\
       \delta_i &= \textrm{Filter width along dimension } i \textrm{ for anisotropic grids}
   \end{aligned}

The anisotropic derivative is used to define the eddy viscosity as

.. math::

   \begin{aligned}
       \tau_{ij} &= -2 \nu_t \widetilde{S}_{ij} \\
       \nu_t &= \frac{- (\hat{\partial}_k \widetilde{u}_i) (\hat{\partial}_k \widetilde{u}_j) \widetilde{S}_{ij}}{ (\partial_l \widetilde{u}_m) (\partial_l \widetilde{u}_m) }
   \end{aligned}


AMD model (for ABL)
^^^^^^^^^^^^^^^^^^^

The eddy viscosity is calculated using an anisotropic derivative with a
different filter width in each direction

.. math::

   \begin{aligned}
       \hat{\partial}_i &= \sqrt{C} \delta_i \partial_i \textrm{ for } i=1,2,3 \\
       C &= 1/3 \textrm{ Poincare coefficient for } 2^{nd} \textrm{ order gradient} \\
       \delta_i &= \textrm{Filter width along dimension } i \textrm{ for anisotropic grids}\\
       \beta &= g/\Theta_0 \textrm{ Gravity constant over reference temperature}
   \end{aligned}

The anisotropic derivative is used to define the eddy viscosity as

.. math::

   \begin{aligned}
       \tau_{ij} &= -2 \nu_t \widetilde{S}_{ij} \\
       \nu_t &= \frac{- (\hat{\partial}_k \widetilde{u}_i) (\hat{\partial}_k \widetilde{u}_j) \widetilde{S}_{ij} +  \beta (\hat{\partial}_k \widetilde{w}) (\hat{\partial}_k (\widetilde{\Theta} - \langle {\widetilde{\Theta}} \rangle) )  }{ (\partial_l \widetilde{u}_m) (\partial_l \widetilde{u}_m) } \\
       \tau_{\theta j} &= -2 D_e \frac{\partial \widetilde{\Theta}}{\partial x_j} \\
       D_e &= \frac{- (\hat{\partial}_k \widetilde{u}_i) (\hat{\partial}_k \widetilde{\Theta}) \partial_i \widetilde{\Theta} }{(\partial_l \widetilde{\Theta}) (\partial_l \widetilde{\Theta})}
   \end{aligned}

- **Unit tests**

There is a simple unit test for both :math:`\nu_t` and :math:`D_e` in
``unit_tests/turbulence/test_turbulence_LES.cpp`` under
``test_AMD_setup_calc``.

Non-linear Sub-grid Scale Model 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The non-linear model extends the Smagorinsky model by including an extra term computed from the strain and vorticity rate. 
The modification proposed by `Branco (JFM 1997) <https://doi.org/10.1017/S0022112096004697>`_ and implemented in WRF (`Mirocha et. al (MWR 2010) <https://doi.org/10.1175/2010MWR3286.1>`_) 
is the model considered. The sub-grid scale stress tensor is calculated as follows: 

 .. math::
    M_{ij}= -(C_s \Delta)^2 
    [
      2(2S_{mn}S_{mn})^{1/2}S_{ij}+C_1(S_{ik}S_{kj}-\frac{1}{3}S_{mn}S_{mn} \delta_{ij})
      +C_2(S_{ik}R_{kj}-R_{ik}S_{kj})
    ]

Here :math:`S_{ij}` is the strain-rate tensor and :math:`R_{ij}` is the vorticity rate tensor. The model constants are: 
:math:`C_s=[8*(1+C_b)/27\pi^2]^{1/2}`, :math:`C_1=C_2=960^{1/2}C_b/7(1+C_b)S_k`, :math:`S_k=0.5`, and :math:`C_b=0.36`.  

The default length scale of :math:`L=C_s\Delta` causes over-prediction of the mean wind speed profiles. To avoid this over-prediction, the
length scale is modified as proposed by `Senocak et. al (BLM 2007) <https://doi.org/10.1007/s10546-007-9181-x>`_

.. math::
   L=(1-\exp(-z/H))^2(\frac{\kappa z}{\phi_M})^2+(\exp(-z/H))^2(C_s \Delta)^2

Here the term :math:`H=1.5 dz` specifies the location at which the length scale switches to :math:`L=C_s\Delta` and :math:`\phi_M`
is the atmospheric stability function. Currently, the implementation for the stability function uses a single global value. 
The implementation of the non-linear model is split into two parts. The subgrid-scale viscosity term is directly used 
within the Kynema-SGF diffusion framework. The last two terms in :math:`M_{ij}` are added as source-terms in the momentum equation. 


