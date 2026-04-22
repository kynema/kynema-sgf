.. _wall_models:

Wall models
-----------
The wall models described in this section are implemented in Kynema-SGF for
running wall-bounded flows.

Monin-Obukhov Similarity Theory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Monin-Obukhov similarity theory is used for wall boundary conditions for ABL simulations. The exact
calculation of :math:`\tau_{i3}` in the horizontal directions depends on the SGS model used, but the following calculations for the friction velocity :math:`u_\tau` and surface heat flux `q` are common across the models.

.. math::
    u_\tau = \frac{\kappa \overline{s}}{\ln \left(\frac{z_b}{z_0}\right) - \psi_m}
    
where :math:`s` is the horizontal wind speed :math:`s = \sqrt{u_{1}^2+ u_{2}^2}`, :math:`\theta_w`
is the wall temperature, :math:`\kappa` is the von Karman constant, and :math:`z_0` is the surface roughness length and :math:`z_b` is the reference height (default is the first cell center). The
:math:`\overline{\phantom{l}.\phantom{l}}` operator indicates a horizontal plane
average.  The quantities :math:`\psi_m, \psi_h` are computed using the Monin-Obukhov similarity law
following the calculations in `ven der Lann et al <https://doi.org/10.1002/we.2017>`_ and `Dyer (1974)` formulation  for unstable stratification (:math:`z_b/L < 0`):

.. math::
    \begin{align}
        \psi_m &= 2\ln \left(\frac{1+x}{2}\right) + \ln \left(\frac{1+x^2}{2}\right) - 2 \arctan{x} + \frac{\pi}{2}, x = \left(1 - \beta_m\frac{z_b}{L}\right)^{\frac{1}{4}} \\
        \psi_h &= \ln \left( \frac{1 + y}{2}\right), y = \left(1 - \beta_h \frac{z_b}{L}\right)^{\frac{1}{2}},
    \end{align}

and for stable stratification (:math:`z_b/L > 0` ):

.. math::
    \begin{align}
        \psi_m &= -\gamma_m \frac{z_b}{L},\\
        \psi_h &= -\gamma_h \frac{z_b}{L},
    \end{align}

where :math:`L = -\frac{u_\tau^3 \theta_0}{\kappa g q}` is the Monin-Obukhov length and :math:`\beta_m, \beta_h, \gamma_m, \gamma_h` are model constants. Kynema-SGF uses :math:`\beta_m = \beta_h = 16` and :math:`\gamma_m = \gamma_h = 5`.

Log-law wall model
~~~~~~~~~~~~~~~~~~

This wall model computes the local :math:`u_\tau` from the velocity at
the first grid cell, and uses this to compute the shear stress, which is
then used as a boundary condition.

The log law:

.. math::
    :label: loglaw

    u_{\mathrm{mag}} = u_\tau \left(\frac{1}{\kappa}\log\left(\frac{u_\tau z}{\nu}\right) + B\right). 

Given a horizontal velocity magnitude
:math:`u_{\mathrm{mag}} = \sqrt{u^2 + v^2}` at
:math:`z = z_{\mathrm{ref}}`, :math:`u_\tau` can be computed using a
non-linear solve to satisfy equation :eq:`loglaw`.

In Kynema-SGF Newton-Raphson iterations are used with a convergence
criterion of :math:`\lvert u_\tau^{n+1} - u_\tau^n \rvert < 10^{-5}`.
For this, derivative of
:math:`\frac{\partial u_{\mathrm{mag}}}{\partial {u_\tau}}` is used,

.. math:: \frac{\partial u_{\mathrm{mag}}}{\partial {u_\tau}} = \left(\frac{1}{\kappa}\left(1+\log\left(\frac{u_\tau z_{\mathrm{ref}}}{\nu}\right)\right) + B\right)

.. math:: u_\tau^{n+1} = u_\tau^{n} - \left(u_\tau^n \left(\frac{1}{\kappa}\log\left(\frac{u_\tau^n z_{\mathrm{ref}}}{\nu}\right) + B\right) - u_{\mathrm{mag}}\right)/\frac{\partial u_{\mathrm{mag}}}{\partial {u_\tau}}.

Finally, the shear stress is calculated as,

.. math::

   \begin{aligned}
       \tau_{xz} &= u_\tau^2 \frac{u}{u_\mathrm{mag}} \\
       \tau_{yz} &= u_\tau^2 \frac{v}{u_\mathrm{mag}}
   \end{aligned}

Constant stress model
~~~~~~~~~~~~~~~~~~~~~

NOTE: This wall model will be ill-posed unless combined with a Dirichlet
boundary condition on the other wall, :math:`\langle u \rangle` can
drift by a constant otherwise.

This is a trivial wall model, where the shear stresses are specified as
constants. For a pressure gradient driven channel,

.. math::

   \begin{aligned}
       u_\tau^2 &= -\frac{\mathrm{d} P}{\mathrm{d} x} \\
       \tau_{xz} &= u_\tau^2 \\
       \tau_{yz} &= 0
   \end{aligned}

Schumann model
~~~~~~~~~~~~~~

NOTE: This wall model will be ill-posed unless combined with a Dirichlet
boundary condition on the other wall, :math:`\langle u \rangle` can
drift by a constant otherwise.

This model is a modified version of the constant stress model, where the
fluctuations from a reference height :math:`z_\mathrm{ref}` are used to
add fluctuations in the shear stress.

.. math::

   \begin{aligned}
       u_\tau^2 &= -\frac{\mathrm{d} P}{\mathrm{d} x} \\
       \tau_{xz} &= u_\tau^2 \frac{u}{\langle u_\mathrm{mag} \rangle} \\
       \tau_{yz} &= u_\tau^2 \frac{v}{\langle u_\mathrm{mag} \rangle}
   \end{aligned}

where, :math:`\langle u_\mathrm{mag} \rangle` is the planar average of
:math:`u_{\mathrm{mag}} = \sqrt{u^2 + v^2}` at :math:`z_\mathrm{ref}`.

Symmetric wall boundary
~~~~~~~~~~~~~~~~~~~~~~~

This is a boundary condition to for flows with a symmetry across the
z direction (example: *half-channel* simulations) at the centerline.

.. math::

   \begin{aligned}
       \tau_{xz} &= 0 \\
       \tau_{yz} &= 0 \\
       w &= 0
   \end{aligned}

Dynamic wall model (Wave model)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This wall model is used to calculate the stress due to moving surfaces,
like ocean waves. It aims to introduce wave phase-resolving physics 
at a cost similar to using the Log-law wall model, without the need of using
wave adapting computational grids. The model was developed by `Ayala et al. (2024) <https://doi.org/10.1007/s10546-024-00884-8>`_.

.. math:: \tau_{i3} = \frac{1}{\pi}|(\boldsymbol{u-C}) \cdot \boldsymbol{\hat{n}}|^2|\boldsymbol{\nabla} \eta|^2 \, \hat{n}_i  \, \text{H} \Bigl[ (u_j-C_j)\frac{\partial \eta}{\partial x_j} \Bigr] \, + \, \tau^{visc}_{i3}, \quad i = 1,2.

The first component gives the form drag due to ocean waves, where :math:`\boldsymbol{C}`
is the wave velocity vector, :math:`\eta` is the surface height distribution and
:math:`\hat{\boldsymbol n} = \boldsymbol{\nabla} \eta /|\boldsymbol{\nabla} \eta|`. The
second component (:math:`\tau^{visc}_{i3}`) is the stress due to unresolved effects,
like viscous effects. For this component, the ``Log-law wall model`` is used.
