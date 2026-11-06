Forest Model
--------------
The forest model provides an option to include the drag from forested regions to be included in the momentum equation. The 
drag force is calculated as follows: 

.. math::

   F_i= - C_d L(x,y,z) U_i | U_i |


Here :math:`C_d` is the coefficient of drag for the forested region and :math:`L(x,y,z)` is the leaf area density (LAD) for the 
forested region. A three-dimensional model for the LAD is usually unavailable and is also cumbersome to use if there are thousands
of trees. Two different models are available as an alternative: 

.. math::
   L=\frac{LAI}{h}

.. math:: 
   L(z)=L_m \left(\frac{h - z_m}{h - z}\right)^n  exp\left[n \left(1 -\frac{h - z_m}{h - z}\right )\right]

Here :math:`LAI` is the leaf area index and is available from measurements, :math:`h` is the height of the tree, :math:`z_m` is the location 
of the maximum LAD, :math:`L_m` is the maximum value of LAD at :math:`z_m` and :math:`n` is a model constant with values  6 (below :math:`z_m`) and 0.5 
(above :math:`z_m`), respectively. :math:`L_m` is computed by integrating the following equation: 

.. math::
   LAI = \int_{0}^{h} L(z) dz 

The simplified model with uniform LAD is recommended for forested regions with no knowledge of the individual trees. LAI values can be used from 
climate model look-up tables for different regions around the world if no local remote sensing data is available. 

In addition to the analytical LAD profiles above, Kynema-SGF also supports a point-cloud representation of each forest. In this
approach, each forest is described by a separate file containing scattered samples of leaf area density,

.. math::

   \left(x_p, y_p, z_p, L_p\right), \qquad p = 1, \ldots, N_f,

where :math:`N_f` is the number of samples for a given forest and :math:`L_p` is the LAD value associated with sample point
:math:`p`. Each forest file is paired with its own drag coefficient :math:`C_d`.

The point-cloud implementation reconstructs a local LAD field directly from the scattered samples instead of assuming a prescribed
vertical profile. To avoid applying forest drag in empty regions between disconnected samples, the sample points are first projected
onto the :math:`x-y` plane and a convex hull is formed. LAD interpolation is only performed for cell centers whose
:math:`(x,y)` location lies inside this convex hull.

For a cell center located at :math:`\boldsymbol{x}=(x,y,z)`, the :math:`k` nearest point-cloud samples from the same forest are selected,
where :math:`k` is controlled by the input parameter ``ForestDrag.point_neighbors``. Let :math:`d_p` denote the Euclidean distance
from the cell center to a selected sample point,

.. math::

   d_p = \sqrt{(x-x_p)^2 + (y-y_p)^2 + (z-z_p)^2}.

If the cell center coincides with a sample point to within the regularization tolerance
``ForestDrag.point_interp_eps`` :math:`= \varepsilon`, then the LAD is taken directly from that sample,

.. math::

   L(x,y,z) = L_p.

Otherwise, the LAD is reconstructed with inverse-distance weighting,

.. math::

   L(x,y,z) = \frac{\sum_{p=1}^{k} w_p L_p}{\sum_{p=1}^{k} w_p},

with weights

.. math::

   w_p = \frac{1}{\sqrt{d_p^2 + \varepsilon^2}}.

This regularization avoids singular behavior when a cell center is extremely close to a sample point.

The current implementation also restricts interpolation in the vertical direction using the selected neighboring samples. Let
:math:`z_{\max}^{(k)}` be the maximum height among the :math:`k` nearest selected samples. Then the interpolated LAD is only applied when
the bottom of the cell is below that local canopy height, i.e.,

.. math::

   z - \frac{\Delta z}{2} \le z_{\max}^{(k)}.

If the cell center lies outside the projected convex hull, or if the cell is above the local canopy height inferred from the nearest
samples, the LAD is taken to be zero and no forest drag is applied.

The point-cloud model is useful when remote-sensing products or preprocessed canopy datasets provide spatially varying LAD samples for
individual forest patches. Compared with the simplified uniform or analytical vertical-profile models, this approach allows the drag field
to vary in all three spatial directions while still using a compact set of scattered sample points.

Canopy turbulence
~~~~~~~~~~~~~~~~~

The canopy drag removes mean kinetic energy from the flow. Part of it becomes
turbulence in the wakes of leaves and branches at scales much smaller than the
mesh, which then dissipates quickly (the spectral short cut). With the
one-equation ``KLAxell`` model, ``ForestDrag.canopy_tke = true`` adds these
effects to the ``KransAxell`` TKE source,

.. math::

   S_k = C_d L \left( \beta_p |U|^3 - \beta_d |U| k \right),

where :math:`C_d L` is the ``forest_drag`` field, :math:`\beta_p` is the
fraction of the drag work converted to TKE (``ForestDrag.canopy_beta_p``,
default 1) and :math:`\beta_d` sets the short-circuit dissipation
(``ForestDrag.canopy_beta_d``, default 4, following
:cite:t:`green:1992`, :cite:t:`liu-etal:1996` and :cite:t:`sanz:2003`). The sink is a relaxation of :math:`k` at
the rate :math:`c = \beta_d C_d L |U|`. It is integrated exactly over a time
step,

.. math::

   S_{k,\mathrm{sink}} = -k \, \frac{1 - e^{-c \Delta t}}{\Delta t},

which tends to :math:`-c \, k` for :math:`c \Delta t \ll 1` and never
removes more than the local :math:`k` in one step, so a dense canopy cannot
drive :math:`k` negative.

Inside a dense canopy the eddies are limited by the foliage rather than by
the distance to the ground. ``KLAxell`` computes its length scale
:math:`l = \lambda \kappa z / (\lambda + \kappa z)` from the height above the
terrain, which overestimates the mixing near the canopy top. With
``ForestDrag.canopy_tke = true`` the length scale is limited by the canopy
drag length :math:`L_c = 1 / (C_d L)`,

.. math::

   l \leftarrow \min \left( l, \ \alpha L_c \right)
   \qquad \text{where } C_d L > 0,

where :math:`\alpha` is ``ForestDrag.canopy_length_alpha``. :cite:t:`harman-finnigan:2007`
give the mixing length in the canopy as
:math:`l = 2 \beta^3 L_c` with :math:`\beta = u_* / U_h` the ratio of the
friction velocity to the wind speed at the canopy top, so
:math:`\alpha = 2 \beta^3`. The default :math:`\alpha = 0.04`
(:math:`\beta \approx 0.27`) was chosen from a sweep of
:math:`\alpha` = 0.02 to 0.07 in a one-dimensional analogue of the
:cite:t:`shaw-schumann:1992` large-eddy simulation (canopy depth ratio, drag
coefficient and leaf area density profiles of their LAI = 2 and 5 cases):
it gave the smallest combined error in the mean wind and momentum flux
profiles for both canopies, and the runs settle at
:math:`\beta = u_*/U_h` of 0.25 to 0.27, consistent with
:math:`2 \beta^3 \approx 0.04`. The typical :math:`\beta = 0.3` would give
0.054, which leaves the wind in a sparse canopy too strong. The eddy viscosity, the shear production and the
buoyancy production of ``KLAxell`` are proportional to :math:`l` and are
rescaled with it; the dissipation :math:`C_\mu^3 k^{3/2} / l` of the
``KransAxell`` source uses the limited :math:`l`. Cells without forest drag
are not changed.

The canopy terms are off by default and are independent of the momentum drag:
with ``ForestDrag.canopy_tke = false`` the forest only acts through
``ForestForcing``, as before.

Roughness representation
~~~~~~~~~~~~~~~~~~~~~~~~

When the canopy is not resolved by the mesh, ``ForestDrag.model = roughness``
replaces the canopy drag by a roughness length in the footprint of each forest.
The terrain wall model then applies the log law in the first cell above the
ground,

.. math::

   u_* = \frac{\kappa \, U_1}{\ln \left( z_1 / z_0 \right) - \psi_m},

with the forest :math:`z_0` given directly or as a fraction of the tree
height, :math:`z_0 = c \, h`. A typical value is :math:`c \approx 0.1`. The
flow sees the forest as a rough surface at the ground: no displacement height
:math:`d` is applied, so the mean profile above a tall canopy is the log law
:math:`\ln(z/z_0)` rather than :math:`\ln((z-d)/z_0)`.

