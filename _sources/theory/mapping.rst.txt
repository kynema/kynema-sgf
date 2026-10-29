.. _mapping:

Mapping definition
==================

Let the mapping of stretched mesh to uniform mesh be defined by :math:`{\mathbf X} = \Phi (\Xi) = (x(\chi), y(\eta), z(\xi))`. For demonstration purposes, we assume mapping only in the :math:`z`-direction. Note that
the normals to faces continue to be aligned with the grid axes.

We define the velocity

.. math:: \overline{U} = T U

where

.. math:: T = J [\nabla_\Xi \Phi]^{-1}

and

.. math:: J = {\rm det} [\nabla_\Xi \Phi]

In our case :math:`T` is

.. math::

   \begin{bmatrix} 
   z_\xi & 0 & 0 \\
   0 & z_\xi & 0\\
   0 & 0 & 1 \\
   \end{bmatrix}, 
   \hspace{0.5in} J = ( z_\xi),
   \hspace{0.5in} \overline{U}  = (z_\xi u, \; z_\xi  v, \; w).

We can write the mesh spacing in :math:`{\mathbf X}`-space as

.. math:: ( x_\chi dx, y_\eta dy, z_\xi dz) =  ( dx, dy, z_\xi dz)

where :math:`(dx,dy,dz)` are the mesh spacings in
:math:`{\mathbf \Xi}`-space.

Handy identities
----------------

.. math:: \nabla_X = \frac{1}{J} T^T \nabla_\Xi

.. math:: \nabla_X \cdot A = \frac{1}{J} \nabla_\Xi \cdot (T A) = \frac{1}{J} \nabla_\Xi \cdot (\overline{A})

for any vector field :math:`A`

Governing equations in stretched mesh space
-------------------------------------------

.. math::

   JU_t + (\overline{U^{MAC}}\cdot \nabla_\Xi) \; U = 
    J(F_b -\frac{1}{\rho} \frac{T^{T}}{J} \nabla_\Xi \; p) + \frac{\mu}{\rho} \nabla_\Xi \cdot (\frac{1}{J} T T^T \nabla_\Xi U) ,

.. math:: J \rho_t + \nabla_\Xi \cdot (\rho \overline{U^{MAC}}) = 0,

.. math:: J (\rho s)_t +  \nabla_\Xi \cdot (\rho s \overline{U^{MAC}}) = S,

.. math:: \nabla_\Xi  \cdot \overline{U}= 0   \;\;\;,

MAC Projection
--------------

In the MAC projection we want to solve

.. math:: \nabla_X  \cdot (\frac{1}{\rho} \; \nabla_X p) = \nabla_X \cdot U^{pred}

then set

.. math:: U^{MAC} = U^{pred} - \frac{1}{\rho} \nabla_X \; p

Mapping the above to uniform coordinates yields

.. math:: \nabla_\Xi  \cdot (\frac{1}{J \rho} T T^T \nabla_\Xi \; p) = \nabla_\Xi \cdot \overline{U^{pred}}

The MAC projection from AMReX will then return

.. math:: \overline{U^{MAC}} = \overline{U^{pred}} - \frac{1}{\rho} T (\frac{1}{J} T^T \nabla_\Xi p) = \overline{U^{pred}} - \frac{1}{J \rho} T T^T \nabla_\Xi \; p

Nodal Projection
----------------

Analogously to the MAC projection we want to solve

.. math:: \nabla_\Xi  \cdot (\frac{1}{J \rho} T T^T \nabla_\Xi \; p) = T^T \nabla_\Xi \cdot ( U^{n+1,*} + \frac{1}{J} T^T \nabla_\Xi p^{n-1/2} ) = \nabla_\Xi \overline{U}^{n+1,*}

where here :math:`\overline{U}^{n+1,*} = T U^{n+1,*}` is cell-centered
with
:math:`U^{n+1,*} = U^{n+1,*} + \frac{1}{J} T^T \nabla_\Xi p^{n-1/2}`.
Note, although :math:`p` is defined on the nodes, :math:`\nabla p` is
cell-centered.

We construct the divergence and gradient operators with a finite element
approach, so the mapping enters in here via integrals over cell volumes.
Thus we will only need :math:`z_\xi` at cell centers; this is the same
metric term used in the MAC projector.

As with the MAC projection we want to create
:math:`\overline{U}^{n+1,*} = T U^{n+1,*}`
as the velocity field we will send in to the nodal projector.

The nodal projection will return

.. math:: \overline{U}^{n+1} = \overline{U}^{n+1,*} - \frac{1}{J \rho} T T^T \nabla_\Xi \; p

Note that – unlike in the MAC projection – we want :math:`U^{n+1}`
rather than :math:`\overline{U}^{n+1},` so after the projection we must
define

.. math:: U = T^{-1} \overline{U}

Initial conditions & computing the timestep
-------------------------------------------

The initial conditions for any problem as well as any other global computations
such as evaluating the CFL or the adaptive timestep size should be done so in the
stretched coordinates. 

Exceptions to the current implementation
----------------------------------------

- As a first pass, the stretched mesh capability is implemented only for the MOL scheme.

- The stretched mesh capability has been tested and verified only for single-level AMR meshes.

- The mesh stretching capability requires the use of multi-component velocity solves by setting ``velocity_diffusion.use_tensor_operator = false``

- Efforts are underway to extend the capability beyond laminar physics.
