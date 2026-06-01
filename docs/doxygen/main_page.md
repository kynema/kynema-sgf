# Kynema-SGF API documentation {#mainpage}

This document is intended for developers who want to understand the C++ code
structure and modify the codebase and, therefore, assumes that the reader is
familiar with the installation, compilation, and execution steps. If you are new to
Kynema-SGF and haven't installed/used Kynema-SGF previously, we recommend starting
with the [user manual](https://kynema.github.io/kynema-sgf/user/user.html) that provides a detailed
overview of the installation process as well as general usage.

## How to use this API guide?

This section provides a brief overview of the organization of this API guide so
as to enable readers to quickly find the sections that they are interested in.
The source code documentation is organized in [sections](modules.html) that
divides the codebase into logical groups. We recommend that you start
[here](modules.html) and navigate to the sections that you are interested in.

Kynema-SGF is built on top of the [AMReX
library](https://amrex-codes.github.io/amrex/). The \ref core "core data structures" 
provide higher-level abstractions on top of AMReX data structures.
We recommend reading the [AMReX basics
chapter](https://amrex-codes.github.io/amrex/docs_html/Basics.html) to
familiarize yourself with the core AMReX terminology and concepts. Once you have
read that chapter, read the \ref core "Kynema-SGF core"
documentation and familiarize yourself with the concept of \ref
kynema_sgf::Field "Field" and \ref kynema_sgf::FieldRepo "FieldRepo" (see \ref
fields) in Kynema-SGF as these are used quite heavily everywhere in the code. Two
other global data structures that are used frequently are \ref kynema_sgf::CFDSim
"CFDSim" and \ref kynema_sgf::SimTime "SimTime". `CFDSim` represents the
simulation environment and holds references to the mesh, the field repository,
time instance, the registered \ref physics "physics" instances, the \ref eqsys
"equation systems", and \ref utilities "post-processing and I/O" utilities.
`SimTime` holds all attributes related to time within the code and determines
when to advance the simulation, exit, or write outputs.

### Source code organization

Upon successful download/clone, the base repository (`kynema-sgf`) has source code
is organized in subdirectories described below:

- `src` -- C++ source files. All code is located within this directory
- `unit_tests` -- Unit-tests for individual modules/classes
- `cmake` -- Functions, utilities used during CMake configuration phase
- `docs` -- User manual (Sphinx-based) and Doxygen files
- `submods` -- Third-party libraries and dependencies
- `test` -- Regression tests and associated infrastructure
- `tools` -- Miscellaneous post-processing scripts and other utilities

When developing new features, we strongly recommend creating a unit-test and
develop features incrementally and testing as you add capabilities. Unit-tests
are also a good way to explore the usage of individual components of the code.

## Contributing

Kynema-SGF is an open-source code and we welcome contributions from the community.
Please consult the [developer
documentation](https://kynema.github.io/kynema-sgf/developer/index.html) section
of the user manual to learn about the process of submitting code enhancements,
bug-fixes, documentation updates, etc.

## License

Kynema-SGF is licensed under BSD 3-clause license. Please see the
[LICENSE](https://github.com/Kynema/kynema-sgf/blob/development/LICENSE) included in
the source code repository for more details.

