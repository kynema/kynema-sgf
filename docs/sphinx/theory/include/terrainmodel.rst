.. _terrainmodel:

Terrain Model
--------------

An immersed boundary forcing method (IBFM) is used to represent the terrain. In this method,
the effect of the terrain is modeled using a forcing term in the momentum and energy equation.
Two implementations are available: the original ``TerrainDrag`` physics with a binary blanking
of cells, and the ``ImmersedTerrain`` physics with a partial terrain fraction per cell. Both
follow the immersed body force method of
`Muñoz‐Esparza, Domingo, et al. (JAMS 2020) <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020MS002141>`_,
which prescribes only the forcing inside the body; the wall functions, stability corrections and
time integration described below are additions made in this code.

Binary blanking (``TerrainDrag``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The forcing term in the momentum equation is given by:

.. math::

   F_i = - \beta C_d u_i | u_i |

Here :math:`\beta` is the volume fraction of the cell covered by terrain, :math:`C_d` is a drag
term  and :math:`u_i` is the wind speed. In ``TerrainDrag`` the volume fraction is
computed as a 0 or 1 using a simple nearest cell algorithm at each grid level, which turns slopes
into a staircase. The calculation of the drag coefficient term and the forcing term for the energy
equation can be found in the reference above.

The original formulation is designed for low Reynolds number cases and does not include a
method for applying a wall function. We propose the use of a forcing function to include
the wall effects.

First, compute the friction velocity from the cell above the terrain-adjacent cell, whose
center is :math:`1.5\,\Delta z` above the wall. The velocity of the terrain-adjacent cell itself
is not used because that is the cell being forced:

.. math::

   u_*= |u_h[k+1]| \frac {\kappa}{\log [1.5 \Delta z/z_0] - \psi_m(1.5 \Delta z / L)}

The expected wind speed at cell k, whose center is :math:`0.5\,\Delta z` above the wall, is

.. math::

   |u_n|= \frac{u_*}{\kappa} \left[ \log (0.5 \Delta z/z_0) - \psi_m(0.5 \Delta z / L) \right]

with :math:`\psi_m` the Monin-Obukhov stability function for a single prescribed Obukhov
length :math:`L` (neutral when it is not specified). The forcing term is computed as

.. math::

   F_i= - \frac {|u[k]| \hat{c} - |u_n|\hat{l}} {\tau}

Here :math:`\hat{c}=(1,1,1)` is the existing normal vector from the grid and :math:`\hat{l}=(ux,uy,0)/|u_n|` is the value
from the log law.

.. image:: ./images/terrain_normal.png
   :align: center
   :width: 30%

Partial terrain fraction (``ImmersedTerrain``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``ImmersedTerrain`` replaces the binary blanking with the fraction :math:`\beta \in [0,1]` of each
cell occupied by terrain (the fraction of the cell column below the terrain height). It also
stores an integer mask (0 fluid, 1 solid, 2 surface) and the terrain height. Surface cells are the
partially filled cells and the fluid cells sharing a face with a solid cell on any of the six
sides, so that the vertical faces of steep terrain and buildings are treated as walls as well.
The mask can be used directly by ``FieldRefinement`` to refine the surface band.

**Immersed drag.** The drag inside the terrain is applied by ``ImmersedDragForcing`` as a linear
relaxation toward zero velocity at the rate :math:`C = w C_d/\Delta z`, with the drag weight
:math:`w = \beta` (``drag_weight = fraction``) or :math:`w = 1` where the cell center is inside the
terrain and 0 elsewhere (``drag_weight = center``). Because the term is
linear in :math:`u`, it is integrated exactly over the time step,

.. math::

   \frac{\partial u_i}{\partial t} = - C_\mathrm{eff}\, u_i, \qquad
   C_\mathrm{eff} = \frac{1 - e^{-C \Delta t}}{\Delta t} \le \min\left(C, \frac{1}{\Delta t}\right),

which is stable and monotone for any :math:`C \Delta t` and removes the need for drag limiters.
The explicit form :math:`-C u_i` overshoots at cold start whenever :math:`C \Delta t > 2`, which
happens routinely on the first CFL-limited step.

**Implicit projection.** With the drag applied explicitly, the velocity inside the body is
reset to zero by the source term and then re-created by the pressure projection, which treats the
body as fluid and adds :math:`\Delta t \nabla p / \rho` back every step. Under a CFL time step this
residual is proportional to :math:`\Delta x`, so the approach to zero inside the body is first order
regardless of how the interface is represented. The option ``ImmersedTerrain.implicit_projection``
applies the drag implicitly together with the pressure,

.. math::

   u^{n+1} = \frac{u^{**} - \Delta t \nabla p / \rho}{1 + \beta C \Delta t},
   \qquad \nabla \cdot u^{n+1} = 0
   \;\Rightarrow\;
   \nabla \cdot \left( \frac{\Delta t}{\rho (1 + \beta C \Delta t)} \nabla p \right)
   = \nabla \cdot \frac{u^{**}}{1 + \beta C \Delta t},

that is, the terrain behaves as a fluid of density :math:`\rho (1 + \beta C \Delta t)` in the
nodal and MAC projections. The residual inside the body becomes independent of the time step and
equal to :math:`\nabla p / (\rho C)`; increasing :math:`C_d` in proportion to :math:`1/\Delta z`
then gives second-order convergence. On the laminar immersed box case the fitted orders of the
mean speed inside the body were 1.3 (binary), 1.1 (partial fraction, explicit drag), 1.0
(implicit projection, fixed :math:`C_d`) and 1.9 (implicit projection, :math:`C_d \propto 1/\Delta z`).

**Implicit drag in the diffusion solve.** With ``implicit_projection`` the terrain cells must
also be held at rest inside the implicit diffusion solve: the operator uses the same effective
density :math:`\rho(1 + \beta C \Delta t)` as the coefficient of the time-derivative term (with the
right-hand side kept at the plain density), so that momentum diffusing into the terrain is damped
there rather than accumulated and removed by the projection. Without this the terrain cells float
to a fraction of the fluid velocity during each solve, the interface flux is reduced, the effective
wall sits about one cell too low, and the coupled scheme becomes unstable when
:math:`\nu \Delta t / \Delta z^2` exceeds about 10.

**Diffusive flux at the interface.** The diffusion operator evaluates the viscous flux at a
fluid/solid face as :math:`\mu(u_k - 0)/\Delta_f`, since the interior is at rest. This is a wall
stress with the wrong length scale: the flux to a no-slip wall at distance :math:`d_1` is
:math:`\mu u_k/d_1`. ``ImmersedTerrain`` therefore multiplies the coefficient of every wall face
by :math:`\Delta_f/d_1`, which reproduces the no-slip flux at the true wall position. With
``center`` weighting the wall
face is the one next to a solid cell and :math:`d_1 = z_k - h` on the bottom face. With
``fraction`` weighting the partial cell is the wall cell: its velocity stands for the centroid of
its fluid part, so the bottom face uses :math:`d_1 = (1-\beta)\Delta z/2`, the face above it is
scaled by :math:`2/(2-\beta)` for the longer distance to the next cell center, only a cell entirely
inside the terrain makes a face a wall, and the partial cell takes no drag.

**Laminar channel verification.** For plane Poiseuille flow with an immersed flat bottom wall
between cell centers and a no-slip top wall, ``implicit_projection`` with the no-slip wall flux
converges to the exact parabolic profile at second order with both drag weights (L2 error falling by a factor of four per refinement, independent of the
sub-cell wall position), provided the drag coefficient is large enough that the residual velocity
of the terrain cells, of order :math:`1/(1 + C \Delta t)`, stays below the discretization error, or
is increased with resolution.
