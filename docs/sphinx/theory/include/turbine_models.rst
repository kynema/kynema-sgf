.. _turbine_models:

Turbine models
--------------

The modeling approaches for turbines, namely actuator disk model (ADM)
and actuator line model (ALM), are presented in this section. Although
Kynema-SGF offers some stand-alone ADM implementations for turbines and
ALM pathways for individual wings, the primary methods for turbine
simulations rely on coupling with the modular, open-source,
multifidelity wind-energy tool
OpenFAST :cite:p:`jonkman2018full,openfast-ref` for blade
motion and turbine behavior while Kynema-SGF calculates aerodynamic
interactions.

Actuator disk model
~~~~~~~~~~~~~~~~~~~

The ADM is used to represent wind turbines as forces in the momentum
equation :cite:p:`Madsen1997,Jimenez_2007,Calaf2010`,
spreading the induced forces over the swept area of the rotor (the
disk). In this approach, the overall thrust force, :math:`F_t`, is
generally calculated using the velocity averaged, :math:`\overline{u}`,
over the rotor-swept area, :math:`\frac{\pi}{4} D^2`, and an empirically
derived constant of thrust, :math:`C_t`, as shown by the formula

.. math:: F_t = -\frac{1}{2}\rho \, C_t \overline{u}^2 \frac{\pi}{4} D^2.

The rotation of the turbine also introduces circulation in the flow via
an azimuthal force, which has different representations depending on the
model variant.

The forces are computed at discrete points on the disk surface using the
actuator point :math:`i`, and these discrete forces are then spread
over the turbine disk via a spreading function, :math:`\Gamma_i`, often
a Gaussian :cite:p:`martinez2015`:

.. math::

   \Gamma_i(x_j) = \frac{1}{\epsilon^3\pi^{3/2}} \exp\bigg(-\frac{h_{i,j}^2}{\epsilon^2}\bigg),
       \label{eq:uniform_gauss}

where :math:`x_j` is an arbitrary location in space, :math:`\epsilon` is
the characteristic width of the Gaussian kernel, and :math:`h_{i,j}` is
the distance between the actuator point :math:`i` and the
:math:`x_j` position in space. A more general form of
`[eq:uniform_gauss] <#eq:uniform_gauss>`__ can be rewritten as a tensor
product of contributions from each component of an orthogonal coordinate
system:

.. math::

   \Gamma_i(x_j) = \eta (h_{i,1}) e_1 \otimes \psi (h_{i,2}) e_3 \otimes \vartheta  (h_{i,3} ) e_3,
       \label{eq:nonuniform_gauss}

where :math:`\eta`, :math:`\psi`, and :math:`\vartheta` represent
different spreading kernels along the three principal coordinates
(:math:`e_1`, :math:`e_2`, and :math:`e_3`) of an orthogonal coordinate
system.

Kynema-SGF offers the traditional Gaussian along with a set of projections
constructed from products of normalized linear basis functions,
:math:`\xi`, similar to finite-elements, where the principal directions
of the spreading functions align with the cylindrical coordinates of the
disk. The geometry of the disk can be resolved with far fewer points
because the linear basis functions form a partition of unity. The linear
basis function can be expressed as

.. math::

   \xi_i(\tilde{x},\Delta h) = \frac{\max \left(0, 1-\frac{\vert \tilde{x}_i-\tilde{x} \vert}{\Delta h}\right)}{\Delta h},
       \label{eq:linear-basis}

where :math:`\Delta h` is the distance between the actuator points in
the associated principal direction and :math:`\tilde{x}` is the
principal direction’s abscissa.

Employing the linear basis function, Kynema-SGF offers a spreading
function that is uniform in the :math:`\theta` direction:

.. math:: \Gamma_i(x_j) = \frac{\xi_r(r, \Delta r)}{2 \pi r_i} \:  \frac{1}{\epsilon \sqrt{\pi}} \exp^{-\frac{\left(\tilde{z}_i-\tilde{z}_j\right)^2}{\epsilon^2}},

where :math:`r` is the radial location along the disk and
:math:`\tilde{z}` is in the normal direction. Similarly, another option
is available that discretizes in the :math:`\theta` direction using arc
length as the characteristic length scale:

.. math:: \Gamma_i(x_j) = \xi_r(r, \Delta r) \: \xi_\theta(r\theta, r_i\Delta \theta) \:  \frac{1}{\epsilon \sqrt{\pi}} \exp^{-\frac{\left(\tilde{z}_i-\tilde{z}_j\right)^2}{\epsilon^2}}.

Figure `1 <#fig:adm-disk-projections>`__ illustrates the different ADM
spreading-function approaches.

.. figure:: ./images/adm_disk_projection.png
   :alt: adm-disk-projections
   :width: 66.0%

   Illustration of actuator disk forcing using traditional point
   Gaussian functions (left) and using products of linear basis functions in the
   disk plane (right): (a) and (b) show the footprint of a single
   actuator point; (c) and (d) show the sum of all the points. The red
   dots indicate the location of the actuator points and the color bar
   shows the strength of a normalized force.

In conjunction with these spreading function options, Kynema-SGF includes
three implementations of the ADM, and in each case the wind speed is
sampled at discrete points, spaced evenly across the disk. To accurately
estimate the freestream velocity, sample points are also placed upstream
of the turbine disk. At each sample point, the sampled wind velocity is
obtained via linear interpolation. After the body forces are calculated
according to the chosen method, they are applied to the flow at the
discrete disk points. This means that the number of points for the wind
speed sampling and the modeled body force can be independently
configured for optimal sampling based on the chosen model and scenario.

The simplest of the ADM options in Kynema-SGF is the uniform :math:`C_t`
model, which employs coefficients of thrust specified by the user at
different wind speeds. To represent the loading on the blades, which
directly translates to the force of the blades on the flow, the uniform
:math:`C_t` model applies the forcing term uniformly over each radial
section of the disk. Further complexity is included in the Joukowsky
disk model, described by :cite:`SORENSEN20202259` and
updated in :cite:`Sorenson2023`. This model relies on an
assumption of constant circulation over the disk with corrections for
hub and tip effects, providing analytical expressions for axial and
azimuthal force distributions. The Joukowsky disk approach has been
validated for different wind turbines in diverse operating regimes and,
as an analytical model, offers notable computational efficiency.
Finally, Kynema-SGF features an ADM approach coupled to OpenFAST
:cite:p:`Cheung_2023`, whereby OpenFAST modules use the
sampled wind speeds to calculate wind turbine dynamics, including the
aerodynamics of the turbine blades. After calculating the resulting body
forces, these are distributed back to the fluid domain within Kynema-SGF
by applying an isotropic smoothing kernel.

In some scenarios, the first two approaches may be more expedient, and
they have the advantage of being fully contained within Kynema-SGF.
However, by coupling to OpenFAST, additional aspects of turbine behavior
are included in the modeling framework, such as structural dynamics, and
controllers can be easily incorporated to govern turbine behavior during
the course of a simulation.

Actuator line model
~~~~~~~~~~~~~~~~~~~

The ALM is used to represent wings and individual wind turbine blades as
body forces in the momentum equation
:cite:p:`sorensen2002numerical,troldborg2012comparison,martinez2015`,
dividing the forces along a blade or wing into forces at actuator points
along a line. In this approach, aerodynamic forces such as lift and drag
are computed at each actuator point using a sampled velocity from the
fluid solver and lookup tables for lift and drag coefficients,
:math:`C_l` and :math:`C_d`:

.. math::

   F_l = \frac{1}{2} \rho \, c \, w \, C_l (\alpha) \, U_{ rel}^2,
       ~~~~~~~~~~~~
       F_d = \frac{1}{2} \rho \, c \, w \, C_d (\alpha) \, U_{ rel}^2,

where :math:`\alpha` is the angle of attack, :math:`c` is the local
chord, :math:`w` is the width of the blade element, and :math:`U_{rel}`
is the magnitude of the relative velocity
:cite:p:`Sorensen:2002,martinez2015`. The relative velocity
vector, which contributes to :math:`\alpha` and :math:`U_{rel}`, is
calculated
:math:`\boldsymbol{u}_{rel} = \boldsymbol{u}(\boldsymbol{x}_{i}) - \dot{\boldsymbol{x}}_{i}`,
where :math:`\dot{\boldsymbol{x}}_i` is the velocity of the
actuator point :math:`i`, :math:`\boldsymbol{x}_i` is the position
of the actuator point :math:`i`, and
:math:`\boldsymbol{u}(\boldsymbol{x}_i)` is the sampled fluid velocity
at that point. The lift force and drag force are then added in vector
form and distributed around surrounding grid cells by using a finite
Gaussian kernel:

.. math::

   \boldsymbol{S}_i(\boldsymbol{x}_j) = \frac{\boldsymbol{f}_i} {\epsilon_c \, \epsilon_t \, \epsilon_s \, \pi^{3/2}}
       \exp_{fin} \left(
           - \frac{x_c^2}{\epsilon_c}
           - \frac{x_t^2}{\epsilon_t}
           - \frac{x_s^2}{\epsilon_s}
             \right),
       \label{eq:lift-coord-spreading}

where :math:`\boldsymbol{f}_i=F_{l,i}\hat{e}_l + F_{d,i}\hat{e}_d` is
the aerodynamic force vector of an actuator point :math:`i`,
:math:`x_c`, :math:`x_t`, and :math:`x_s` are the distances from the
center of the actuator point :math:`i` to the point
:math:`\boldsymbol{x}_j` in the chord, thickness, and spanwise
directions, and :math:`\epsilon_c`, :math:`\epsilon_t`, and
:math:`\epsilon_s` are the thickness of the Gaussian in the
corresponding directions :cite:p:`churchfield2017`. A finite
Gaussian kernel is used to avoid calculating unreasonably small force
values over the entire domain, and it is defined using

.. math::

   \exp_{fin}(-h^2) = \left\{\begin{matrix}
   \exp (-h^2) & \textrm{if } h < 4 \\
   0 & \textrm{otherwise}
   \end{matrix}\right. .

The Godunov-based time discretization of Kynema-SGF specifies that forcing
terms, such as :math:`\boldsymbol{S}_i` should be calculated at the half
time step, :math:`n+1/2`. Calculating :math:`\boldsymbol{S}_i` consists
of two parts: the point-force vector :math:`\boldsymbol{f}_i` and the
Gaussian kernel centered at :math:`(x_c,x_t,x_s)`. To get an accurate
value of :math:`\boldsymbol{f}_i`, the velocity field must be sampled at
the location of the actuator point, enabling the computation of the
relative velocity :math:`U_{rel}`. Because the velocity field at
:math:`n+1/2` is not available when the point force needs to be
calculated, we use the latest available information for this step; we
sample the velocity field :math:`\boldsymbol{u}^n` at the location of
the actuator point at time step :math:`n`, which is
:math:`(x_c,x_t,x_s)^n`. Next, the Gaussian kernel translates the point
force to a force field, centered at the actuator point. Because of the
time discretization, we place the Gaussian kernel at the actuator
location at time step :math:`n+1/2`, written as
:math:`(x_c,x_t,x_s)^{n+1/2}`. This location is straightforward to
obtain through the turbine model interface: after sampling the velocity
field and computing the point force, the turbine model advances to
:math:`n+1`, providing the actuator point locations at :math:`n+1`.
Actuator point locations at :math:`n+1/2` are calculated through a
simple arithmetic average:

.. math:: (x_c,x_t,x_s)^{n+1/2} = \frac{1}{2}\left(x_c^n+x_c^{n+1},x_t^n+x_t^{n+1},x_s^n+x_s^{n+1} \right). \label{eq:alm_loc}

Because the accuracy of ALM relies on sampling the velocity at the
location of the bound vortex :cite:p:`Martinez2017`, and this
vortex is created from the inclusion of the source term
:math:`\boldsymbol{S}_i`, the relationship between velocity sampling and
force placement influences the accuracy of the method. Accuracy issues,
such as dependence of the actuator force on the time step size and
incorrect power estimation for the turbine, arise when
`[eq:alm_loc] <#eq:alm_loc>`__ is not used for the location of the
actuator force, neglecting the time-staggered nature of Kynema-SGF.

To simulate wind turbines using ALM, Kynema-SGF relies on OpenFAST for the
turbine representation. In this coupled framework, Kynema-SGF supplies
sampled velocities at actuator points to OpenFAST at the beginning of
every time step, and OpenFAST returns updated actuator point forces and
locations after progressing the turbine model forward in time. These
forces and locations are incorporated into the Kynema-SGF numerical
algorithm through the actuator force implementation. By using OpenFAST,
including aero-elastic deformations and changes in wind turbine operation
is straightforward. These effects are computed by OpenFAST modules as
the simulation evolves and are passed on to Kynema-SGF as changes to the
actuator point locations.

Kynema-SGF offers a test bed for advanced actuator lines, providing
variations to force distribution and point placement. The advanced ALM
features include anisotropic Gaussian body forces to better represent
wind turbine blades :cite:p:`churchfield2012large` and the
filtered lifting line correction to obtain consistent results across
different grid resolutions and :math:`\epsilon` values
:cite:p:`martinez-tossas_meneveau_2019, dag2020, martinez2023`.
