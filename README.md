# Kynema-SGF 

[Documentation](https://kynema.github.io/amr-wind) | [Nightly test dashboard](http://my.cdash.org/index.php?project=Exawind) 

[![Powered by AMReX](https://amrex-codes.github.io/badges/powered%20by-AMReX-red.svg)](https://amrex-codes.github.io/amrex/) [![Build Status](https://github.com/Exawind/amr-wind/workflows/Kynema-SGF-CI/badge.svg)](https://github.com/Exawind/amr-wind/actions) [![OpenSSF Best Practices](https://www.bestpractices.dev/projects/9284/badge)](https://www.bestpractices.dev/projects/9284)


Kynema-SGF is a massively parallel, block-structured adaptive-mesh, incompressible
flow solver for wind turbine and wind farm simulations. The codebase is a
wind-focused fork of [incflo](https://github.com/AMReX-Codes/incflo). The solver
is built on top of the [AMReX library](https://amrex-codes.github.io/amrex).
AMReX library provides the mesh data structures, mesh adaptivity, as well as the
linear solvers used for solving the governing equations. Kynema-SGF is actively
developed and maintained by a dedicated multi-institutional team from [Lawrence
Berkeley National Laboratory](https://www.lbl.gov/), [National Renewable Energy
Laboratory](https://nrel.gov), and [Sandia National
Laboratories](https://sandia.gov).

The primary applications for Kynema-SGF are: performing large-eddy simulations
(LES) of atmospheric boundary layer (ABL) flows, simulating wind farm
turbine-wake interactions using actuator disk or actuator line models for
turbines, and as a background solver when coupled with a near-body solver (e.g.,
[Nalu-Wind](https://github.com/kynema/nalu-wind)) with overset methodology to
perform blade-resolved simulations of multiple wind turbines within a wind farm.
For offshore applications, the ability to model the air-sea interaction effects
and its impact on the ABL characteristics is another focus for the code
development effort. As with other codes in the
[Exawind](https://github.com/kynema) ecosystem, AMR-wind shares the following
objectives:

- an open, well-documented implementation of the state-of-the-art computational
  models for modeling wind farm flow physics at various fidelities that are
  backed by a comprehensive verification and validation (V&V) process;

- be capable of performing the highest-fidelity simulations of flow fields within
  wind farms; and 

- be able to leverage the high-performance leadership class computing
  facilities available at DOE national laboratories.

## Part of the WETO Stack

Kynema-SGF is primarily developed with the support of the U.S. Department of Energy and is part of the [WETO Software Stack](https://NatLabRockies.github.io/WETOStack). For more information and other integrated modeling software, see:
- [Portfolio Overview](https://NatLabRockies.github.io/WETOStack/portfolio_analysis/overview.html)
- [Entry Guide](https://NatLabRockies.github.io/WETOStack/_static/entry_guide/index.html)
- [High-Fidelity Modeling Workshop](https://NatLabRockies.github.io/WETOStack/workshops/user_workshops_2024.html#high-fidelity-modeling)

## References
To cite Kynema-SGF and to learn more about the methodology, use the following [journal article](https://doi.org/10.1002/we.70010) as well as the [ExaWind reference](https://doi.org/10.1002/we.2886):

```
@article{amrwind2025,
    author = {Kuhn, Michael B. and {Henry de Frahan}, Marc T. and Mohan, Prakash and Deskos, Georgios and Churchfield, Matthew and Cheung, Lawrence and Sharma, Ashesh and Almgren, Ann and Ananthan, Shreyas and Brazell, Michael J. and {Martinez-Tossas} Luis A. and Thedin, Regis and Rood, Jon and Sakievich, Philip and Vijayakumar, Ganesh and Zhang, Weiqun and Sprague, Michael A.},
    title = {Kynema-SGF: A performance-portable, high-fidelity flow solver for wind farm simulations},
    journal = {Wind Energy},
    volume = {28},
    number = {5},
    pages = {-},
    doi = {10.1002/we.70010},
    url = {https://onlinelibrary.wiley.com/doi/full/10.1002/we.70010},
    eprint = {},
    year = {2025}
}

@article{exawind2024,
    author = {Sharma, Ashesh and Brazell, Michael J. and Vijayakumar, Ganesh and Ananthan, Shreyas and Cheung, Lawrence and deVelder, Nathaniel and {Henry de Frahan}, Marc T. and Matula, Neil and Mullowney, Paul and Rood, Jon and Sakievich, Philip and Almgren, Ann and Crozier, Paul S. and Sprague, Michael},
    title = {ExaWind: Open-source CFD for hybrid-RANS/LES geometry-resolved wind turbine simulations in atmospheric flows},
    journal = {Wind Energy},
    volume = {27},
    number = {3},
    pages = {225-257},
    doi = {10.1002/we.2886},
    url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/we.2886},
    eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/we.2886},
    year = {2024}
}
```

## Documentation

Documentation is available at https://kynema.github.io/amr-wind, which 
includes a walkthrough tutorial, a user manual, notes on theory,
and tips for developers. We also provide a developer-focused API
documentation at the same link. You can either
browse the docs online by following the links, or you can generate them locally
after downloading the code. Please follow the instructions in user manual to
build documentation locally.

## Compilation and usage

Kynema-SGF is built upon the [AMReX library](https://amrex-codes.github.io/amrex).
A snapshot of the AMReX library is distributed along with the Kynema-SGF source
code as a `git-submodule`. In addition to the AMReX library, you will require a
modern C++ compiler that supports the C++17 standard. Users wishing to execute
the code on high-performance computing (HPC) systems will also need MPI
libraries installed on their system. The code can also be compiled using MPI+X, 
where X can be OpenMP for CPU shared memory parallelism,
CUDA to target NVIDIA GPUs, ROCM for AMD GPUs, or SyCL for Intel GPUs.

### Contributing, reporting bugs, and requesting help

To report issues or bugs please [create a new
issue](https://github.com/Exawind/amr-wind/issues/new) on GitHub.

We welcome contributions from the community in form of bug fixes, feature
enhancements, documentation updates, etc. All contributions are processed
through pull-requests on GitHub. Please refer to the 
[coding guidelines](https://kynema.github.io/amr-wind/developer/coding_guidelines.html) as
a reference for the best practices currently used to develop Kynema-SGF.

To pass the code format checking, use this command with a new version of `clang-format`:
```
find src unit_tests \( -name "*.cpp" -o -name "*.H" -o -name "*.h" \) -exec clang-format -i {} +
```

Please acknowledge as a publication co-author any developer that has
significantly contributed to implementing or improving specific
capability that was used for that publication.

### User discussion, feedback, and community support

The development team manages a mailing list for Kynema-SGF users. Invites for quarterly user meetings,
along with occasional announcements, are sent to this list.
Quarterly meetings provide development updates and a forum for discussion and feedback.
If you would like to join this mailing list, please send a request to amr-wind-maintainers@groups.nrel.gov,
and we will be happy to add your email address. Our maintainers email is also available for direct
inquiries about Kynema-SGF, but the GitHub page (issues, discussions, pull requests) is preferred
for the majority of questions.

## Versioning and tags

Kynema-SGF uses a type of semantic versioning to help users navigate different versions of the code, 
which are labeled with GitHub tags. These tagged versions are not exhaustive, and they adhere to
the following convention. Given a version number MAJOR.MINOR.PATCH:
1. MAJOR version for changes to input file compatibility for key aspects of the solver, when a key model is changed to significantly affect results of simulations, when a major new capability is added
2. MINOR version for when a significant feature is added (in a backward compatible manner), accumulation of smaller features, or changes to input file compatibility for less central aspects of the solver (e.g., post-processing, forcing terms)
3. PATCH version for backward compatible bug fixes

## License

Kynema-SGF is licensed under BSD 3-clause license. Please see the
[LICENSE](https://github.com/Exawind/amr-wind/blob/development/LICENSE) included in
the source code repository for more details.

