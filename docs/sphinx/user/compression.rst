.. _compression:

Compressing field output files
==============================

When writing out field output files, sometimes called plot files or plt files, it is possible to output the files with HDF5 and compress these with ZFP. To enable this capability, Kynema-SGF must be compiled with ``KYNEMA_SGF_ENABLE_HDF5`` and ``KYNEMA_SGF_ENABLE_HDF5_ZFP``. If using a spack installation of Kynema-SGF, simply add ``+hdf5`` to the spec to enable these compression features. The input file :ref:`documentation <inputs_io>` documents the input options to output these files. The `AMReX HDF5 documentation <https://amrex-codes.github.io/amrex/docs_html/IO.html#hdf5-plotfile-compression>`_ further documents valid options.

While generating the compressed files is straightforward, it takes a little work to read these into traditional visualization software (e.g., Paraview, Visit). However, the HDF5 files are easier to open with Python than traditional AMReX plotfiles.

ParaView
--------

To read ZFP compressed HDF5 files into ParaView, ParaView must have the environment variable ``HDF5_PLUGIN_PATH`` populated and pointing to the HDF5 ZFP plugin location. Within a spack environment with the ``h5z-zfp`` package installed:

.. code-block:: console

    $ export HDF5_PLUGIN_PATH=$(spack location -i h5z-zfp)/plugin
    $ paraview

Then the user can load the file as usual with the ``VisitChomboReader``.

Convert to numpy files
----------------------

To convert the HDF5 files to plain numpy data files, we have provided a script to do so in `tools` directory of this repository. It requires the ``h5py`` and ``hdf5plugin`` Python packages.

.. code-block:: console

    $ python convert_amrex_hdf5_plotfile.py -f plt117501-hdf5-zfp-0.0000001

The resulting file can then be opened in Python with numpy.

yt
--

To read HDF5 files into `yt <https://yt-project.org>`_, the HDF5 file attributes need to be modified. We have provided a script to do so in `tools` directory of this repository. It requires the ``h5py`` and ``hdf5plugin`` Python packages.

.. code-block:: console

    $ python modify_hdf5_attributes.py -f plt00000 plt00010

The resulting file can then be opened in yt.
