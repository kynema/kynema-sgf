.. _governing_equations:

Governing equations
-------------------


Kynema-SGF solves the LES formulation of the incompressible Navier-Stokes
equations with additional source terms and a transport equation for potential
temperature to model atmospheric boundary layer flows. With
:math:`\widetilde{.}` denoting the spatial filtering operator, the governing
equations in Cartesian coordinates using Einstein notation are

.. math::

   \begin{aligned}
     \frac{\partial \widetilde{u_j}}{\partial x_j} & = 0,\label{eqn:ns-les-cont}\\
     \frac{\partial \widetilde{u_i}}{\partial t} +
     \frac{\partial \widetilde{u_i} \widetilde{u_j}}{\partial x_j} &=
     - \frac{1}{\rho}{\frac{\partial \widetilde{p}}{\partial x_i}}
     - {\frac{\partial \tau_{ij}}{\partial x_j} }
     + \nu \frac{\partial^2 \widetilde{u_i}}{\partial x_j \partial x_j}
     + {C_i}
     + {B_i}
     + {F_{i}},\label{eqn:ns-les-mom}\\
     \frac{\partial \widetilde{\theta}}{\partial t} +
     \frac{\partial \widetilde{u_j} \widetilde{\theta}}{\partial x_j} &= 
     - \frac{\partial \tau_{\theta j}}{\partial x_j} + \frac{\nu}{\mathrm{Pr}} \frac{\partial^2 \widetilde{\theta}}{\partial x_j \partial x_j}, \label{eqn:pot-temp-les}
   \end{aligned}

where :math:`x_i` denotes the coordinate in direction :math:`i`,
:math:`\widetilde{u_i}` is the filtered velocity, :math:`\widetilde{p}` is the
pressure, and :math:`\widetilde{\theta}` is potential temperature; :math:`\nu`
is the molecular viscosity and Pr is the laminar Prandtl number;
:math:`\tau_{ij}` and :math:`\tau_{\theta j}` are the subgrid stress and
heat flux terms, representing the interaction with the unresolved
quantities, defined as

.. math::

   \begin{aligned}
       \tau_{ij} &= \widetilde{u_i u_j} - \widetilde{u_i}\widetilde{u_j} \\
       \tau_{\theta j} &= \widetilde{\theta u_j} - \widetilde{\theta}\widetilde{u_j}.
   \end{aligned}

:math:`C_i` is the contribution of the Coriolis forces due to Earth’s
rotation and is modeled as

.. math:: C_i = -2 \epsilon_{ijk}\Omega_j{u_k},

where :math:`\Omega_j` is the Earth’s angular velocity. :math:`B_i` is
the buoyancy term, modeled using a Boussinesq approximation as

.. math:: B_i = -g_i \beta \left( \theta - \theta_0 \right),

where :math:`\theta_0` is the reference ambient temperature,
:math:`\beta` is thermal expansion coefficient, approximated as
:math:`\beta \approx 1 / \theta_0`, and :math:`g_i` is the gravity
vector. Through this model, the density :math:`\rho` is constant and
uniform in the governing equations. :math:`F_i` represents other force
terms that may be applied in a simulation. This includes geostrophic
forcing, which drives the flow to a desired horizontal mean velocity.
:math:`F_i` also includes the forces introduced by turbines via actuator
methods.

For ease of notation, other parts of this documentation will drop the
explicit use of filtering operator :math:`\widetilde{.}` with the understanding
that the discussion is on the filtered quantities. Additionally, the notation
:math:`(x,y,z) = (x_1, x_2, x_3)` and :math:`(u,v,w) = (u_1, u_2, u_3)`,
will be used to denote the coordinate directions and velocity components in the
corresponding directions.
