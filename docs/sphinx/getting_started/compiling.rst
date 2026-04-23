.. _compiling:

Compiling using Spack, with Exawind-Manager
===========================================

The first step before using Kynema-SGF is compilation of the software on the computing resources
you intend to use. We recommend the use of Spack to streamline this process because of its ability
to handle code dependencies and keep track of hardware variations.

Spack is a "package management tool designed to support multiple versions and configurations of 
software on a wide variety of platforms and environments." To make it easier for users to harness
Spack, we provide the repository `Exawind-Manager <https://github.com/Exawind/kynema-manager>`_.
Exawind-Manager is a custom application of the Spack-Manager tool, which uses Spack environments 
to manage software builds. From the 
`Spack-Manager documentation <https://sandialabs.github.io/spack-manager/user_profiles/developers/developer_spack_minimum.html>`_,
"These environments are similar to 
Conda environments in concept, but they benefit from reusing software that you've built in previous environments. 
As such it is recommended that you maintain a single instance of Spack-Manager to organize and curate your builds, 
and create new environments when you want to start a new development project." 

This walkthrough is meant for common Kynema-SGF workflows, and to avoid verbosity we provide limited 
details about the functionality of Spack and Exawind-Manager. Some additional details are included in 
dropdown sections, but these can be ignored if they are hard to understand.
If you would rather avoid the use of Spack and handle dependencies manually, please refer to
the section of the user manual outlining the process of building from the source code using CMake: :doc:`../user/build`.

.. collapse:: Further (external) installation examples

    Further installation examples are provided in the Spack-Manager documentation, 
    including `Snapshot Developer Workflow Example <https://sandialabs.github.io/spack-manager/user_profiles/developers/snapshot_workflow.html>`_
    and `Spack-Manager abbreviated example <https://sandialabs.github.io/spack-manager/user_profiles/developers/developer_workflow.html#quick-start>`_.
    However, it should be noted that when using Exawind-Manager, the environment variable EXAWIND_MANAGER should be used in place of
    the SPACK_MANAGER variable, pointing to the location of the Exawind-Manager directory.

|

Set up Exawind-Manager
----------------------

To begin the process, clone Exawind-Manager using the following command, which includes its submodules through the ``--recursive`` option.

.. code-block:: console

    export REPO_DEST=<fill in the directory>
    cd ${REPO_DEST}
    git clone --recursive https://github.com/Exawind/kynema-manager.git

Because the kynema-manager directory will house any dependencies that need to be downloaded, it is best to put it in a
location with sufficient memory.

Set up Spack environment
------------------------

Then activate it:

.. code-block:: console

    export EXAWIND_MANAGER=${REPO_DEST}/kynema-manager
    source ${EXAWIND_MANAGER}/start.sh && spack-start

It is helpful to put the above commands as a function in your bash environment so they can be easily 
called in the future when code needs to be altered and recompiled.

Before installing anything, we can obtain info about the codes that Spack supports using ``spack info``:

.. code-block:: console
    
    spack info kynema-sgf

and 

.. code-block:: console
    
    spack info openfast

This command shows the available versions for the specified package as well as different flags (variants)
that can be included, along with the default values for the flags. Turning on a bool-type variant is
represented by a ``+`` symbol, and turning off is denoted with a ``~`` symbol.

For this walkthrough, we will include variants that are commonly relevant to wind turbine simulations
that rely on the actuator-line method (ALM). The next step is to create the environment. 

.. code-block:: console

    mkdir ${REPO_DEST}/env_amrwind_openfast && cd ${REPO_DEST}/env_amrwind_openfast
    quick-create-dev -d . -s kynema-sgf@main+openfast+netcdf%compiler openfast@3.5.3+rosco%compiler

This ``quick-create-dev`` command has flags selected so that that Kynema-SGF will work with OpenFAST,
Kynema-SGF can save certain files using NetCDF, and OpenFAST will compile with
the turbine controller package ROSCO. By executing this command, the environment is set up and activated,
and the Kynema-SGF and OpenFAST repositories are cloned to the environment directory. 

.. collapse:: Details on quick-create-dev

    Environments can be associated with directories (using the ``-d`` *option) or with a name (using the*
    ``-n`` option). Environments associated with directories tend to be more navigable for development 
    because named environments create and use directories within the Exawind-Manager directory.

    The repositories cloned by Exawind-Manager are shallow clones, and do not automatically have any commit 
    history. If you would like to compile an older version of a code using a different commit, you can retrieve
    the commit history using the command ``git fetch --unshallow`` within the repository and then check out 
    any past commit that you may need. After choosing a different commit, be sure to run ``git submodule update``
    to modify the submodules to correspond to the chosen commit.
    
    If you do not want
    to clone new copies of Kynema-SGF or OpenFAST and instead want to use other, already-cloned repositories:
    after making the environment directory, create symbolic links to the cloned repositories in the environment
    directory, ensuring that the name of the links match the name of the repository. Then create the
    environment with ``quick-create-dev``. If you use your own cloned repositories, be aware that this 
    approach puts the version-based dependency checks into your own hands, though.
    
    The fact we specified
    main and develop branches when we created the environment does not mean that the code in these repositories 
    must be on the main and develop branches, respectively. These references communicate to Spack a grouping of 
    dependencies for each code. In many cases, using different commits or even your own fork for ExaWind codes will 
    not change their dependencies, and so the specification of main or master is typically the correct spec for 
    whatever version of the code you are using. OpenFAST compatibility can vary more from version to version, though.


The choice of compiler depends on the machine you are using. On Kestrel, the compiler for CPUs is ``oneapi``,
and on a Mac, the recommended compiler is ``apple-clang``. Omitting ``%compiler`` allows Spack to choose
the default compiler based on the machine being used, but this should not be omitted on Kestrel due to
its GPUs and CPUs both being associated with the same machine identification.

The next step is to have Spack compile everything. If you are using a machine that has dedicated compute nodes,
now is the time to get an interactive job on a compute node. Once on the node, the Spack environment will
need to be activated again. If you are using a machine that does not have dedicated compute resources for compiling,
e.g., a laptop, you can begin a new terminal shell to continue along using the following instructions.
Alternatively, the command ``despacktivate`` will deactivate the Spack environment to put your current terminal in a similar state.

Compile and load
----------------

To activate the Spack environment now, first activate Spack by repeating these commands from before:

.. code-block:: console

    export EXAWIND_MANAGER=${REPO_DEST}/kynema-manager
    source ${EXAWIND_MANAGER}/start.sh && spack-start

Now, the Spack environment can be activated.

.. code-block:: console

    cd ${REPO_DEST}/env_amrwind_openfast
    quick-activate .

Finally, the Spack compilation is

.. code-block:: console

    spack install

This step takes a long time and generates a lot of output text. Spack first determines which packages are required by dependencies 
and then compiles and installs them, downloading when necessary. If this command is interrupted, it can be resumed by following the
same process of activating the environment and repeating the install command.

.. collapse:: Details on spack install and modifying environments

    The process of determining which packages are required is known as concretizing. When ``spack install`` is called in an
    environment that has not been concretized, the concretize step is automatically included. However, these can be called separately,
    as ``spack concretize`` will perform this part of the installation by itself. If the specifications (specs) of an environment
    are changed after concretization, this step may need to be forced to overwrite the preexisting environment using ``spack concretize -f``.
    Environments can be modified by editing the spack.yaml file or by using the ``spack rm <spec>`` command to remove specs (e.g., 
    ``kynema-sgf@main``) and ``spack add <spec>`` to add specs (e.g., ``kynema-sgf@main+hypre``).

After the code is compiled, the executables can be located within build-spack directories inside the package directories, and each
package build has its own hash. Instead of referencing these locations directly to use the executables, Spack provides a command
to add them to the path, enabling the executable to be used directly. When the spack environment is active, use

.. code-block:: console

    spack load kynema-sgf

to make executables from Kynema-SGF directly available. To verify that the package was loaded correctly, type

.. code-block:: console

    spack find --loaded

which will display all the loaded packages.

.. _rosco-dyn-lib:

Find ROSCO dynamic library
--------------------------

On using ROSCO: OpenFAST requires the location of the ROSCO library file (either ``libdiscon.so`` (Linux) 
or ``libdiscon.dylib`` (Mac)) as an argument within the ServoDyn input file.
During the spack install command, the location of the installed 
packages are printed to the screen. After installation is complete, these can be listed again more briefly
by repeating the spack install command. To find the location of the ROSCO library, look for "rosco" among
the listed locations. If ``<spack opt path>/rosco-<hash>`` is the directory provided by spack install,
the ``libdiscon.so`` or ``libdiscon.dylib`` file will be located within ``<spack opt path>/rosco-<hash>/lib/``.

The path to this library file will come into play when setting up the turbine simulation.

|

Go to the next step: :doc:`../walkthrough/precursor`
