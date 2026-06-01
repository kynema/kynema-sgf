Solving physics on a stretched AMReX mesh
------------------------------------------

Often for simulations involving walls, (e.g., channel flows, complex terrains etc.) it is desirable to have finer mesh near the wall which gradually coarsens only in the wall-normal direction. Consequently, modifications within Kynema-SGF are underway to support solving a non-uniformly spaced mesh while still using most of AMReX's machinery directed at uniformly-spaced Cartesian meshes. The governing equations solved on a non-uniform stretched mesh are further explained below -

.. toctree::
   :glob:
   :maxdepth: 2

   mapping.rst
