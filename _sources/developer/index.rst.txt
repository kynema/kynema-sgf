.. _dev-docs:

Developer Documentation
=======================

This section is meant for potential developers who are interested in modifying
and extending the Kynema-SGF source code for their own use cases.

Kynema-SGF is built on top of `AMReX library
<https://amrex-codes.github.io/amrex/docs_html/>`_. Users are strongly
recommended to read through the AMReX documentation and understand the basic
AMReX concepts before jumping into the Kynema-SGF source code.

The `Basics section
<https://amrex-codes.github.io/amrex/docs_html/Basics_Chapter.html>`_ provides a
thorough overview of the basic data structures and ways to interact with these
structures. The `GPU section
<https://amrex-codes.github.io/amrex/docs_html/GPU_Chapter.html>`_ provides an
overview of the AMReX GPU strategy and the higher-level functions (e.g.,
``parallel-for`` abstractions) available to write GPU-ready code within
Kynema-SGF. The `Linear Solvers section
<https://amrex-codes.github.io/amrex/docs_html/LinearSolvers_Chapter.html>`_
provides an overview of the multi-level multigrid (MLMG) solvers used to solve
the various linear systems within Kynema-SGF.

.. toctree::
   :maxdepth: 3

   documentation
   ../doxygen/html/index
   unit_testing
   regression_testing
   verification
   coding_guidelines
