.. _flather_boundary:

Flather open boundary condition
-------------------------------

The Flather boundary condition :cite:p:`flather:1976` is a radiation-type
open boundary for free-surface
flows. It allows surface gravity waves generated inside the domain to leave through
a lateral boundary while still driving the flow toward an externally specified
state. Compared to a simple extrapolation (Neumann) outflow, it avoids the spurious
reflection of long waves back into the domain, and compared to a pure Dirichlet
inflow, it does not over-constrain the interior solution. Comparisons of open
boundary conditions for regional tidal simulations have found this formulation to
be among the better performing choices :cite:p:`carter-merrifield:2007`.

Kynema-SGF uses the condition in a depth-integrated form, which makes it
well suited to the volume-of-fluid representation of the free surface described in
:ref:`multiphase`. All quantities used by the boundary condition are column
integrals taken along the vertical direction, so the boundary condition considers
the transport of the liquid phase rather than focusing on individual cells in isolation.
Though the Flather condition is formulated in two dimensions, it must be applied in the
three-dimensional domain of Kynema-SGF, which involves scaling the local velocity 
with depth-integrated quantities, enabling the preservation of the vertical velocity profile.
This is established practice for regional oceanic models, in which the
Flather condition is applied to the depth-integrated mode while the vertical
structure of the flow is treated separately :cite:p:`marchesiello:2001`.

Depth-integrated quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a given column of cells normal to a lateral boundary, the liquid height and the
depth-integrated normal velocity (a volumetric flux per unit width) are

.. math::

   h = \sum_k \alpha_k\, \Delta z, \qquad
   (uh) = \sum_k u_k\, \alpha_k\, \Delta z,

where :math:`\alpha_k` is the volume fraction of liquid and :math:`u_k` is the
boundary-normal velocity component in cell :math:`k`. Because the volume fraction
appears as a weight, cells that contain only gas make no contribution. These sums
are accumulated across all levels of the AMR hierarchy so that the resulting
profiles are independent of the grid decomposition and of the local refinement.

Two interior column integrals are accumulated separately: :math:`(uh)^{\rm int}_{\rm liq}`
over cells that are entirely liquid, and :math:`(uh)^{\rm int}_{\rm mix}` over
interfacial cells that contain both phases. The corresponding integral in the
boundary (ghost) cells, :math:`(uh)^{\rm ext}`, represents the externally
prescribed state, and :math:`h^{\rm ext}` is its liquid height.

Boundary formulation
~~~~~~~~~~~~~~~~~~~~

The target depth-integrated flux at the boundary is

.. math::

   (uh)^{\rm target} = (uh)^{\rm ext} \pm c \left( h^{\rm int} - h^{\rm ext} \right),
   \qquad c = \sqrt{g\, h^{\rm int}},

where :math:`c` is the shallow-water wave speed, :math:`g` is the magnitude of the
vertical component of :input_param:`incflo.gravity`, and the sign is negative on a
low-side boundary and positive on a high-side boundary. The difference in liquid
height between the interior and the boundary is what radiates the outgoing wave.

Rather than imposing a uniform velocity, the interior velocity profile is rescaled
so that its depth integral matches the target. Only fully liquid cells are
rescaled; interfacial cells are left unchanged, because scaling interfacial cells
have the undesirable consequence of accelerating the gas phase. The scaling factor is

.. math::

   f = \frac{(uh)^{\rm target} - (uh)^{\rm int}_{\rm mix}}
                 {(uh)^{\rm int}_{\rm liq}}.

The interior velocity is also clipped so that it never drives inflow at an
outflow boundary, which keeps the applied profile consistent with the column
integrals used to construct it.

Limiting and edge cases
~~~~~~~~~~~~~~~~~~~~~~~

The scaling above is undefined or poorly conditioned in several situations, so the
implementation falls back to the externally prescribed profile or to a plain
outflow condition in the following cases.

* If the interior column contains no fully liquid cells, there is nothing to
  rescale, and the externally specified profile is used instead.
* If :math:`f` falls outside the range set by
  :input_param:`Flather.min_velocity_scale_factor` and
  :input_param:`Flather.max_velocity_scale_factor`, the externally specified
  profile is used instead. This limit is most relevant during startup, when the
  interior and exterior states can differ substantially and an unbounded scaling
  would produce rapid, nonphysical acceleration.
* If the depth-integrated boundary velocity is essentially zero, the scaling based
  on external quantities is undefined, and the interior profile is retained.
* If the boundary column contains no liquid, the computed wave speed and target
  flux are unreliable, and a standard extrapolation outflow is applied.

Finally, when the depth-integrated boundary flux indicates inflow, the externally
prescribed profile is used. Inflow is assessed from the column integral rather than
from individual cells, so a single cell does not switch the character of the
boundary.
