.. _discretization:

Discretization
--------------

The numerical methodology used to solve the partial differential
equations (PDEs) within Kynema-SGF is documented in `Almgren et
al. (JCP 1998)
<https://ccse.lbl.gov/Publications/almgren/abchw.pdf>`_. Kynema-SGF uses
`AMReX-Hydro
<https://amrex-fluids.github.io/amrex-hydro/docs_html/Schemes.html>`_
for many advection routines. The reader is referred to their
documentation for implementation details.

Time Step -- Godunov
~~~~~~~~~~~~~~~~~~~~

-  Define :math:`U^{MAC,n+1/2}`, the MAC velocity which is used for advection. This velocity is interpolated to the cell faces and extrapolated forward in time using a Taylor expansion. Then it is projected to form a divergence-free velocity field.

.. math:: u_f^{n+1/2} = u^n \pm \frac{\Delta x}{2}\frac{\partial u}{\partial x} &-
          \frac{\Delta t}{2}\left(u^n \frac{\partial u}{\partial x}
          + v^n \frac{\partial u}{\partial y} + w^n \frac{\partial u}{\partial z}\right)
          \\ &+
          \frac{\Delta t}{2}\left(g + \frac{1}{\rho^n}\left(
          -\frac{\partial p^{n-1/2}}{\partial x} + \mu\nabla^2u^n\right) \right)

.. math:: \nabla \cdot \frac{\nabla \phi}{\rho^n} = \nabla \cdot \boldsymbol{u}_f^{n+1/2}

.. math:: \boldsymbol{U}^{MAC,n+1/2} = \boldsymbol{u}_f^{n+1/2} - \frac{\nabla \phi}{\rho^n}

-  Time discretization of momentum governing equation

  - The time step goes from :math:`n` to :math:`n+1` with the right-hand-side at :math:`n+1/2`.
  - The time discretization of :math:`\boldsymbol{\tau}` depends on the method chosen to compute diffusion.

.. math:: (\rho \boldsymbol{U})^{n+1} = (\rho \boldsymbol{U})^n &-
           \Delta t \left( \nabla \cdot (\rho \boldsymbol{U} \otimes \boldsymbol{U}^{MAC})^{n+1/2}
           + \nabla {p}^{n+1/2} \right) \\ &+
           \Delta t \left( \nabla \cdot \boldsymbol{\tau} + \rho^{n+1/2} g \right)

-  Partition the discretized equation into two steps, the Predictor and Applying the Projection. The predicted state, :math:`\ast`, is an approximation to the new-time state, :math:`n+1`.

.. math:: (\rho \boldsymbol{U})^{\ast} = (\rho \boldsymbol{U})^n &-
           \Delta t \left( \nabla \cdot (\rho \boldsymbol{U} \otimes \boldsymbol{U}^{MAC})^{n+1/2}
           + \nabla {p}^{n-1/2} \right) \\ &+
           \Delta t \left( \nabla \cdot \boldsymbol{\tau} + \rho^{n+1/2} g \right)

.. math:: (\rho\boldsymbol{U})^{n+1} = (\rho\boldsymbol{U})^{\ast} + \Delta t \nabla p^{n-1/2} - \Delta t \nabla p^{n+1/2}

- In order to calculate :math:`\rho \boldsymbol{U}^{n+1/2}` within the advection term, the momentum must be interpolated to the faces and extrapolated forward in time a half step, similar to the face velocities involved in the MAC projection. To do this for the momentum, the routine uses the momentum, pressure gradients, source terms, and diffusion terms at :math:`n`, as well as :math:`\boldsymbol{U}^{MAC,n+1/2}`.

- In the case of variable density single-phase or multiphase simulations, the density at :math:`n+1` is found using separate scalar equations, which are solved during the predictor step. Because density has no projection step,

.. math:: \rho^{n+1} = \rho^{\ast}.

Therefore, the equation that applies the new pressure gradient becomes

.. math:: \boldsymbol{U}^{n+1} = \boldsymbol{U}^{\ast} + \frac{1}{\rho^{n+1}}\left(\Delta t \nabla p^{n-1/2} - \Delta t \nabla p^{n+1/2}\right)

- The pressure gradient at :math:`n+1/2` is found by solving the projection

.. math:: \nabla \cdot \frac{1}{\rho^{n+1}} \nabla p^{n+1/2} = \nabla \cdot \left( \frac{1}{\Delta t}
          \boldsymbol{U}^{\ast}+ \frac{1}{\rho^{n+1}} \nabla {p}^{n-1/2} \right)


