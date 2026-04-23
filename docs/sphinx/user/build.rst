.. _build:

Compiling Kynema-SGF
==================

Kynema-SGF is primarily written in C++ and requires a modern C++ compiler (that
supports C++14 standard) to build. This section describes the dependencies of
Kynema-SGF as well as the general procedure to compile and execute Kynema-SGF on
your system. The main dependencies are listed below:

#. Operating system -- Kynema-SGF has been tested on Linux and MacOS operating systems.

#. C++ compiler -- Kynema-SGF requires a C++ compiler that supports the C++14
   standard. The code has been tested on GCC v5.x and higher, LLVM Clang v7.x,
   and Intel 2019 compiler suites.

#. MPI -- `OpenMPI <https://www.open-mpi.org/>`_, `MPICH
   <https://www.mpich.org/>`_, or Intel-MPI.

#. `CMake <https://cmake.org/>`_ -- Configure and build the code

#. `Python <https://python.org>`_ A recent version of python

**Optional dependencies**

#. `NVIDIA CUDA <https://developer.nvidia.com/cuda-zone>`_ version 10 or higher
   required to run on GPUs.

#. `Intel OneAPI
   <https://software.intel.com/content/www/us/en/develop/tools/oneapi.html>`_
   for building with SyCL.

#. `AMD GPU
   <https://www.amd.com/en/products/software/rocm.html>`_
   for building with ROCM.

Building from source
--------------------

#. If you are on an HPC system that provides Modules Environment, load the
   necessary compiler, MPI, and CMake modules. If targeting GPUs, load CUDA
   modules. You can also use scripts from `exawind-builder
   <https://exawind.github.io/exawind-builder>`_.

#. Clone a local copy of the git repository

   .. code-block:: console

      git clone --recursive https://github.com/exawind/amr-wind.git

#. Configure and build

   .. code-block:: console

      mkdir amr-wind/build && cd amr-wind/build
      cmake -DKYNEMA_SGF_ENABLE_TESTS:BOOL=ON ../
      make

   Upon successful build, you will end up with an executable :file:`kynema_sgf` in
   your build directory.

#. (Optional) Test your build

   .. code-block:: console

      ctest --output-on-failure

CMake configuration reference
-----------------------------

Architecture options
~~~~~~~~~~~~~~~~~~~~

.. cmakeval:: KYNEMA_SGF_ENABLE_MPI

   Enable MPI support for parallel executions. Default: OFF

.. cmakeval:: KYNEMA_SGF_ENABLE_OPENMP

   Enable OpenMP threading support for CPU builds. It is not recommended to
   combine this with GPU builds. Default: OFF

.. cmakeval:: KYNEMA_SGF_ENABLE_CUDA

   Enable `NVIDIA CUDA GPU <https://developer.nvidia.com/cuda-zone>`_ builds. Default: OFF

.. cmakeval:: KYNEMA_SGF_ENABLE_SYCL

   Enable `Intel OneAPI SyCL <https://software.intel.com/content/www/us/en/develop/tools/oneapi.html>`_ builds. Default: OFF

.. cmakeval:: KYNEMA_SGF_ENABLE_ROCM

   Enable `AMD ROCM GPU <https://www.amd.com/en/products/software/rocm.html>`_ builds. Default: OFF

Dependencies
~~~~~~~~~~~~~

.. cmakeval:: KYNEMA_SGF_ENABLE_MASA

   Enable MASA library. Default: OFF

.. cmakeval:: KYNEMA_SGF_ENABLE_NETCDF

   Enable NetCDF outputs. Default: OFF

.. cmakeval:: KYNEMA_SGF_ENABLE_OPENFAST

   Enable OpenFAST coupling. Default: OFF

.. cmakeval:: KYNEMA_SGF_OPENFAST_VERSION

   OpenFAST version. Default: 3

Other Kynema-SGF specific options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. cmakeval:: KYNEMA_SGF_ENABLE_DOCUMENTATION

   Enable documentation. Default: OFF

.. cmakeval:: KYNEMA_SGF_ENABLE_TESTS

   Enable CTest testing. Default: OFF

.. cmakeval:: KYNEMA_SGF_TEST_WITH_FCOMPARE

   Enable checking test results against gold files using :program:`fcompare`. Default: OFF

.. cmakeval:: KYNEMA_SGF_REFERENCE_GOLDS_DIRECTORY

   Directory for reference gold files when using :program:`fcompare`. Default: ""

.. cmakeval:: KYNEMA_SGF_SAVE_GOLDS

   Enable saving of the gold files. Default: OFF

.. cmakeval:: KYNEMA_SGF_SAVED_GOLDS_DIRECTORY

   Directory for saving gold files when using :program:`fcompare`. Default: ""

.. cmakeval:: KYNEMA_SGF_TEST_WITH_PYTHON

   Enable performing tests that rely on python, primarily consists of checks 
   for post-processing scripts in tools directory. Default: OFF

General CMake options
~~~~~~~~~~~~~~~~~~~~~~

.. cmakeval:: CMAKE_INSTALL_PREFIX

   The directory where the compiled executables and libraries as well as headers
   are installed. For example, passing
   ``-DCMAKE_INSTALL_PREFIX=${HOME}/software`` will install the executables in
   ``${HOME}/software/bin`` when the user executes the ``make install`` command.

.. cmakeval:: CMAKE_BUILD_TYPE

   Controls the optimization levels for compilation. This variable can take the
   following values:

     ===============  =======================
     Value            Typical flags
     ===============  =======================
     RELEASE          ``-O2 -DNDEBUG``
     DEBUG            ``-g``
     RelWithDebInfo   ``-O2 -g``
     ===============  =======================

   Example: ``-DCMAKE_BUILD_TYPE:STRING=RELEASE``

.. cmakeval:: CMAKE_CXX_COMPILER

   Set the C++ compiler used for compiling the code.

   For Intel SyCL builds (see :cmakeval:`KYNEMA_SGF_ENABLE_SYCL`) this should be
   set to `icpx`.

.. cmakeval:: CMAKE_C_COMPILER

   Set the C compiler used for compiling the code

.. cmakeval:: CMAKE_Fortran_COMPILER

   Set the Fortran compiler used for compiling the code

.. cmakeval:: CMAKE_CXX_FLAGS

   Additional flags to be passed to the C++ compiler during compilation.

.. cmakeval:: CMAKE_C_FLAGS

   Additional flags to be passed to the C compiler during compilation.

.. cmakeval:: CMAKE_Fortran_FLAGS

   Additional flags to be passed to the Fortran compiler during compilation.
