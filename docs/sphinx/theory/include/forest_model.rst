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

