.. spelling:word-list::

   godunov

.. _inputs_sampling_freesurface_sampler:

FreeSurfaceSampler
==================

These examples plot the evolution of the liquid-gas interface as a liquid column 
collapses and flows into the rest of the domain, and they show the velocity 
magnitude (in color) as a function of interface position. At each instant plotted,
the shading of the line color represents the progression of time. Below, the scripts
for processing the data in each available format (``native``, ``netcdf``, and ``ascii``)
are shown, and these scripts are also available in the ``tools`` directory.

Note that the data accesses in the native and ascii examples can be applied to 
other sampler types because the format is identical. However, the NetCDF format
can vary depending on the sampler type, and further detail is provided in that example.

The three examples will generate the same plot, shown here:

.. image:: ./plot_sampling_freesurface.pdf
   :width: 300pt

|

Native format: dam break example
--------------------------------
To generate the data required to replicate this example, run the simulation contained in 
test/test_files/dam_break_godunov with the sampling format set to "native" and the sampling
frequency (actually an interval) set to 30.

.. literalinclude:: ../../../tools/sampling_dam_break_godunov_native.py
   :language: python

NetCDF format: dam break example
--------------------------------
To generate the data required to replicate this example, run the simulation contained in 
test/test_files/dam_break_godunov with the sampling format set to "netcdf" (Kynema-SGF must 
also be compiled with NetCDF) and the sampling frequency (actually an interval) set to 30.

Because the FreeSurfaceSampler tracks the interface location, the sample points change position 
with time, and these changing positions are output as the "points" field. However, samplers that 
do not have moving sample locations do not have the "points" field in the NetCDF dataset. For 
static samplers, the locations are provided solely through the variable "coordinates", which represents
the initial positions of the sampler points. Otherwise, the data accesses in this example can be 
applied to other sampler types.

.. literalinclude:: ../../../tools/sampling_dam_break_godunov_netcdf.py
   :language: python

ASCII format: dam break example
--------------------------------
To generate the data required to replicate this example, run the simulation contained in 
test/test_files/dam_break_godunov with the sampling format set to "ascii" and the sampling
frequency (actually an interval) set to 30.

.. literalinclude:: ../../../tools/sampling_dam_break_godunov_ascii.py
   :language: python
