.. _sampling:

Sampling
========

Kynema-SGF provides a sampling utility that measures field
variables at specific locations within the flow. A variety of
sampler types are available, which determine the locations to 
be sampled and, depending on the type, update the sampling 
locations as time progresses. The sampling data can be output
in AMReX particle native format (``native``), NetCDF format 
(``netcdf``), or AMReX particle ASCII format (``ascii``).
The native and ascii formats are identical for every sampler,
but the netcdf format can include additional output variables
depending on the sampler type.

.. toctree::
   :maxdepth: 2

   sampling_freesurfacesampler.rst
