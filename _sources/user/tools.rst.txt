.. spelling:word-list::

   fcompare
   csv

.. _tools:

Tools reference
===============

This section summarizes the functionality of the auxiliary tools included
in the Kynema-SGF ``tools/`` directory. There are two main parts of this folder. Python
scripts, which do not require compilation, are in ``tools/``. C++ programs, which
are compiled as separate executables when Kynema-SGF is compiled, are saved in ``tools/utilities/``.
After compilation, these executables can be found in the build directory, where each has its own 
folder within ``tools/utilities/`` there.

These capabilities are not meant to be exhaustive and are not maintained as actively as the
solver source code. More tools for interacting with Kynema-SGF data (pre- and post-processing)
can be found in `amr-wind-frontend <https://github.com/Kynema/amr-wind-frontend>`_.

Python scripts
--------------

.. input_param:: amrex_particle.py

    Contains helpful routines for manipulating AMReX particle data.

.. input_param:: amrex_plotfile.py

    Contains helpful routines for manipulating AMReX plotfile data.

.. input_param:: amrex_utils.py

    Contains helpful routines for interacting with AMReX data structures.

.. input_param:: calc_inflowoutflow_stats.py

    Tool to process statistics from precursor ABL runs and provide information to populate certain inputs of a subsequent inflow-outflow simulation.

.. input_param:: convert_amrex_hdf5_plotfile.py

    Converts plotfiles written in HDF5 format to plain numpy data files.

.. input_param:: convert_native_sample_to_time_series.py

    Converts sampler data written in native format to a time series written in ASCII format. This is intended for scenarios when there is a single sampler point of interest, which has to be specified by naming the sampler labels and point index.

.. input_param:: convert_native_sampling_to_structured_ascii.py

    Converts sampler data written in native format to files written in ASCII format. For every sampling folder (i.e. every output step), this sampler creates a file for each sampler group, where each file lists the sampled data in order of the points belonging to that sampler.

.. input_param:: example_plotfile_io.py

    Example script for directly interacting with plotfile data.

.. input_param:: fcompare_particles.py

    Tool to compare native AMReX particle data. This has similar capability to the AMReX ``fcompare`` utility, which compares mesh data written to AMReX plotfiles.

.. input_param:: generate_native_boundary_plane.py

    Tool to generate arbitrary temporal and spatially varying boundary conditions via boundary plane files written in native format.

.. input_param:: generate_native_boundary_plane_header.py

    Tool to generate native format boundary plane header files for arbitrary temporal and spatially varying boundary conditions.

.. input_param:: modify_hdf5_attributes.py

    Modifies HDF5 attributes of files in order to be read into yt.

.. input_param:: native_boundary_plane.py

    Contains helpful routines for manipulating native boundary plane data.

.. input_param:: refine_native_boundary_plane.py

    Apply mesh refinement to a boundary plane file written in native format.

.. input_param:: sampling_dam_break_godunov_ascii.py

    Example script for plotting free surface sampler outputs in ASCII format.

.. input_param:: sampling_dam_break_godunov_native.py

    Example script for plotting free surface sampler outputs in native particle format.

.. input_param:: sampling_dam_break_godunov_netcdf.py
    
    Example script for plotting free surface sampler outputs in NetCDF format.


C++ programs (utilities)
------------------------


.. input_param:: CheckpointToCSV

    Converts checkpoint files to CSV format.

.. input_param:: PlotfileToCSV

    Converts checkpoint files to CSV format.

.. input_param:: coarsen-chkpt

    Reads in a checkpoint file and adds a coarser base level to the existing grid.

.. input_param:: compareMultilevelToReference

    Compares plotfiles (similar to fcompare) when the grid refinements do not exactly match between the two.

.. input_param:: refine-chkpt

    Reads in a checkpoint file and refines it by increasing its base resolution.