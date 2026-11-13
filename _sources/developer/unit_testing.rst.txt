.. _dev-unit-tests:

.. default-domain:: cpp

Unit testing
============

Kynema-SGF uses `GoogleTest
<https://github.com/google/googletest/blob/master/googletest/docs/primer.md>`_
to provide unit-testing capabilities for the source code. Unit tests are built
by default when compiling Kynema-SGF using CMake and can be run using the
:program:`kynema_sgf_unit_tests` executable --- see :ref:`build` for details of
building Kynema-SGF with CMake. This section documents the process of running the
unit tests, the unit test scaffolding facilities available within Kynema-SGF unit
testing infrastructure, and a brief overview of the process of creating new unit
tests.

Running unit tests
------------------

To run the unit test suite, simply execute the :program:`kynema_sgf_unit_tests`
executable at the command prompt. This will execute all the unit tests and print
a summary of the success/failure status for each tests. The executable will also
print a summary of the total number of passed, failed, and skipped tests at the
end of execution. Additional arguments can be used to control the behavior of
execution of unit tests, .e.g., to run a single test or a subset of tests.

.. code-block:: console

   # Show command line arguments and brief help
   ./kynema_sgf_unit_tests -h

   # List available tests
   ./kynema_sgf_unit_tests --gtest_list_tests

   # Run a single test
   ./kynema_sgf_unit_tests --gtest_filter="ABLTest.abl_forcing"

   # Run all tests belonging to a test suite
   ./kynema_sgf_unit_tests --gtest_filter="ABLTest.*"

   # Run all tests beginning with keyword
   ./kynema_sgf_unit_tests --gtest_filter="*.*field*"

   # Run all tests except for those in ABLTest suite
   ./kynema_sgf_unit_tests --gtest_filter="-ABLTest.*"

Basic unit test concepts
------------------------

This section provides a quick overview of basic unit testing concepts for
beginners. We recommend reading at least `GoogleTest Primer docs
<https://github.com/google/googletest/blob/master/googletest/docs/primer.md>`_
to get a better understanding of unit test concepts and how to use GoogleTest to
create new unit tests. More advanced use cases can be found in `GoogleTest
advanced docs
<https://github.com/google/googletest/blob/master/googletest/docs/advanced.md>`_.
Unit tests small functions that test a specific aspect of the code. These tests
are organized into *test suites* that group different tests by category. For
example, all tests related to the definition of *fields* and *field repository*
is organized in the test suite ``FieldRepoTest``. To run just this subset of tests use:

.. code-block:: console

   # Run only tests related to field repository
   ./kynema_sgf_unit_tests --gtest_filter="FieldRepoTest.*"

Assertions
~~~~~~~~~~~

Within a test, we check for expected behavior using `*assertions*
<https://github.com/google/googletest/blob/master/googletest/docs/primer.md#assertions>`_
that test for different conditions. For example, if we were writing a function
``square`` that took one argument, a real number, and returned the square of the
input value, we would write a unit test as shown below:

.. code-block:: c++

   //! Return the square of the number
   double square(double x)
   {
       return x * x;
   }

   TEST(TestSquare, test_square)
   {
       EXPECT_EQ(square(4), 16);
       EXPECT_EQ(square(5), 25);

       // Check negative numbers too
       EXPECT_EQ(square(-4), 16);
       EXPECT_EQ(square(-5), 25);

   }

In the above example, ``EXPECT_EQ`` is an assertion provided by GoogleTest that
allows us to check that the two arguments are equal. We use this to test that
the output of the function matches expected values. Unit tests can be used to
check other expected behaviors of the function rather than just its correctness.
For example, consider a function ``square_root`` that
computes the square root of a given real number, but is supposed to throw a
runtime error if it encounters a negative number. Unit tests allow us to test
both cases: 1. the function produces the correct results for positive numbers
and, 2. it throws an error for negative numbers. The following code shows this example

.. code-block:: c++

   //! Return the square root of a number
   double square_root(double x)
   {
       if (x < 0)
           throw std::runtime_error("Square root requires a positive number");

       return std::sqrt(x);
   }

   TEST(TestSquareRoot, test_sqrt)
   {
       EXPECT_NEAR(square_root(4.0), 2.0, 1.0e-12);
       EXPECT_NEAR(square_root(16.0), 4.0, 1.0e-12);

       // Check sqrt(2) to a lower tolerance (1.0e-3)
       EXPECT_NEAR(square_root(2.0), 1.41421356237, 1.0e-3);

       // Check that negative number creates runtime error
       EXPECT_THROW(square_root(-2.0), std::runtime_error);
   }

Tests and Test Fixtures
~~~~~~~~~~~~~~~~~~~~~~~

GoogleTest supports two types of tests: simple tests, and tests that require
fixtures. Simple tests are tests that test standalone functions that require no
additional data structures for its execution -- see `GoogleTest Simple Test docs
<https://github.com/google/googletest/blob/master/googletest/docs/primer.md#simple-tests>`_.
Test fixtures, on the other hand, allow you to group all necessary data in a
custom test class that can be reused for multiple tests. Defining a new simple
test requires the use of ``TEST()`` macro, and defining a new test with a
fixture requires the use of ``TEST_F()`` macro. In Kynema-SGF, we use test
fixtures extensively to perform actions like creating a mesh and generating some
test fields that will be used to perform the tests.

Unit testing in Kynema-SGF
--------------------------

Unit test files are in :file:`kynema-sgf/unit_tests` directory. All unit test code
is written within the ``kynema_sgf_tests`` namespace. This section describes the
scaffolding available to create unit tests and provides a few examples of unit
tests to help developers write new ones.

.. namespace:: kynema_sgf_tests

Unit test scaffolding
~~~~~~~~~~~~~~~~~~~~~

While unit testing simple functions is straightforward, performing unit tests on
CFD applications is more complicated. Most often this will require the creation
of a test mesh, generating field data structures that will be used to set up and
run the test. This is also complicated by the fact that AMReX creates several
global data structures, e.g., Geometry and ParmParse, that must be properly
reset between each test to ensure a clean environment for each test. Kynema-SGF
provides a few classes that provide the necessary scaffolding to quickly setup
and run tests.

Within the :file:`kynema-sgf/unit_tests` directory, the scaffolding utilities
related to testing are in :file:`ks_test_utils` directory. This section provides
a brief overview of the core classes and their purpose.

``pp_utils`` - ParmParse utilities
```````````````````````````````````

Classes written in Kynema-SGF often require user inputs that are generally read in
from input files through `amrex::ParmParse` (see `docs
<https://amrex-codes.github.io/amrex/docs_html/Basics.html#parmparse>`_).
``pp_utils`` are a set of functions that create skeleton input data that are
used by various classes during initialization. For example, it populates the
problem domain, mesh sizes, etc. so that `amrex::AmrMesh` can be
initialized properly.

.. function:: default_mesh_inputs()

   Populates ParmParse data structure with all the necessary inputs to create an
   AMRMesh instance.

.. function:: default_time_inputs()

   Populates ParmParse data structure with necessary inputs for `kynema_sgf::SimTime`.

.. class:: AmrexTestEnv

:class:`AmrexTestEnv` is a subclass of ``::testing::Environment`` that provides
global setup and teardown actions. This classes is registered with the
GoogleTest environment and is responsible for calling ``amrex::Initialize()``
and ``amrex::Finalize()`` at appropriate times. It also customizes the AMReX
setup by changing a few defaults:

- Disables AMReX's default ``signal_handling`` behavior and restores standard
  C++ exception handling. This is necessary for ``EXPECT_THROW`` type assertions
  to function properly.

- Sets the verbosity such that no messages are printed from AMReX library during
  the execution of unit tests.

- Calling ``ParmParse::Finalize()`` immediately after ``amrex::Initialize()`` so
  that each test can begin with a clean parameter environment.

The last two actions can be overridden by the user for specific invocations of
the unit test executable by providing additional command line arguments. For
example, to set the verbosity:

.. code-block:: console

   ./kynema_sgf_unit_tests amrex.verbose=1

And to keep parameters provided in the command line for use with tests:

.. code-block:: console

   ./kynema_sgf_unit_tests utest.keep_parameters=1

In normal development workflow, users will almost never have to interact with
AmrexTestEnv directly.

.. class:: AmrexTest

:class:`AmrexTest` is a test fixture is derived from ``::testing::Test`` and is
the base class for all the other test fixtures used within Kynema-SGF unit testing
infrastructure. It provides setup and teardown actions that call
``ParmParse::Initialize()`` and ``ParmParse::Finalize()`` actions respectively
to create a clean inputs table environment for each test. The setup/teardown
actions are called before and after a ``TEST_F()`` body is executed. This
fixture does not create a mesh or related data structures, and can be used as a
base fixture for tests that do not require any underlying mesh description.

The following example, taken from one of the unit tests in Kynema-SGF, shows usage
of this test fixture:

.. code-block:: c++
   :linenos:

   // All unit tests are created within the `kynema_sgf_tests` namespace
   namespace kynema_sgf_tests {

   // Anonymous namespace to declare utility functions used within this file
   namespace {

   // Helper function to populate ParmParse with all the necessary inputs used by
   // SimTime class.
   void build_simtime_params()
   {
       amrex::ParmParse pp("time");
       pp.add("stop_time", 2.0);
       pp.add("max_step", 10);
       pp.add("fixed_dt", -0.1);
       pp.add("init_shrink", 0.1);
       pp.add("cfl", 0.45);
       pp.add("verbose", -1);
       pp.add("regrid_interval", 3);
       pp.add("plot_interval", 1);
       pp.add("checkpoint_interval", 2);
   }

   } // namespace

   // Create unique namespace for this test fixture. This is useful to group tests
   // related to SimTime object for filtering during execution.
   class SimTimeTest : public AmrexTest {};

   // This is an example of a unit test that tests SimTime behavior with the
   // AmrexTest test fixture
   TEST_F(SimTimeTest, fixed_dt_loop)
   {
       // Call helper function to populate the defaults
       build_simtime_params();
       {
           // Override defaults to switch to fixed timestep
           amrex::ParmParse pp("time");
           pp.add("fixed_dt", 0.2);
       }

       // Create the object that is being tested
       kynema_sgf::SimTime time;
       time.parse_parameters();

       // Perform tests
       int counter = 0;
       int regrid_counter = 0;
       int plot_counter = 0;
       int chkpt_counter = 0;
       while (time.new_timestep()) {
           time.set_current_cfl(2.0);
           ++counter;

           if (time.write_plot_file()) ++plot_counter;
           if (time.write_checkpoint()) ++chkpt_counter;
           if (time.do_regrid()) ++regrid_counter;
       }
       EXPECT_EQ(counter, 10);
       EXPECT_EQ(plot_counter, 10);
       EXPECT_EQ(chkpt_counter, 5);
       EXPECT_EQ(regrid_counter, 3);

       EXPECT_FALSE(time.write_last_checkpoint());
       EXPECT_FALSE(time.write_last_plot_file());
   }

   } // namespace kynema_sgf_tests


.. class:: AmrTestMesh

:class:`AmrTestMesh` is a concrete implementation of `amrex::AmrCore`
that creates an AMR mesh that can be used with unit testing. In addition to
implementing the basic level data creation methods and refinement routines
``ErrorEst``, it also creates an `kynema_sgf::FieldRepo` instance for
creating and manipulating fields from within unit tests. :class:`AmrTestMesh` is
never directly created within unit tests, instead it is created on-demand
through the test fixture :class:`MeshTest` described next.

.. class:: MeshTest

:class:`MeshTest` is the base test fixture for any test that requires a mesh and
associated field data that will be used by the test. In addition to performing
setup/teardown actions described in :class:`AmrexTest`, it also resets the
default `amrex::Geometry` static data so that different tests can run on
different problem domains prescribed by the test fixture.

Almost all unit tests within Kynema-SGF use :class:`MeshTest` as their base test
fixture. In order to allow grouping tests in logical test suites. The following
example shows the basic usage of this test fixture.

.. code-block:: c++
   :linenos:

   /** Example showing the use of MeshTest test fixture in Kynema-SGF unit tests
    *
    */

   // Kynema-SGF unit test namespace
   namespace kynema_sgf_tests {

   // Create a unique name for this test suite (for grouping and filtering)
   class DemoTest : public MeshTest
   {};

   TEST_F(DemoTest, test_demo_meshtest)
   {
       // Before performing any actions the mesh has to be initialized
       initialize_mesh();

       // Now all data structures are ready for use by the test

       // Access the kynema_sgf::CFDSim object
       auto& sim = sim();

       // Access the simulation time object
       auto& time = sim.time();

       // Access the field repository object
       auto& repo = sim.repo();

       // Access the mesh itself
       auto& mesh = sim.mesh();

       //
       // Perform tests with data
       //

       // By default, field repository must be empty
       EXPECT_EQ(repo.num_fields(), 0);
   }

   } // namespace kynema_sgf_tests

The next example shows a more advanced use case where the user can override
defaults before creating the mesh.

.. code-block:: c++
   :linenos:

   TEST_F(DemoTest, test_meshtest_advanced)
   {
       // This test shows a more advanced example to create intermediate data
       // before generating mesh

       // 1. populate the default parameters
       populate_parameters();

       // 2. Change default mesh resolution
       amrex::Vector<int> ncell{{16, 16, 32}};
       pp.addarr("ncell", ncell);

       // 3. Create the mesh instance
       create_mesh_instance();

       // 4. Declare additional fields
       auto& repo = sim().repo();

       repo.declare_cc_field("velocity", 3, 1, 2);
       repo.declare_cc_field("density", 1, 1, 1);
       repo.declare_nd_field("pressure", 1, 1, 1);

       EXPECT_EQ(repo.num_fields, 3);

       // 5. Create the mesh structure
       initialize_mesh();

       // Check that there is at least 1 level
       EXPECT_EQ(repo.num_active_levels(), 1);
   }


Methods defined by :class:`MeshTest`

.. namespace-push:: MeshTest

.. function:: initialize_mesh()

A test must call :func:`initialize_mesh()` before performing any tests that
require a mesh or associated fields. Behind the scenes,
:func:`initialize_mesh()` performs several actions. It calls
:func:`populate_parameters()`, creates a mesh instance, creates levels from
scratch. After a call to this function, the mesh is ready for use.


.. function:: populate_parameters()

   Populate default parameters necessary for creating an AMRMesh and
   `kynema_sgf::SimTime` objects.


.. function:: create_mesh_instance()

   Create a new AMRMesh instance. This doesn't create the level data from
   scratch yet. That is deferred until :func:`initialize_mesh()` is called.

.. namespace-pop::

Example unit tests
~~~~~~~~~~~~~~~~~~

Following are a list of unit tests available within Kynema-SGF repository that can
be used as starting points for users to write new tests:

`test_simtime.cpp <https://github.com/Kynema/kynema-sgf/blob/development/unit_tests/core/test_simtime.cpp>`_

   Simple unit test example that tests the behavior of
   `kynema_sgf::SimTime`. This test only relies on :class:`AmrexTest` and
   does not require a mesh.

`test_abl_init.cpp <https://github.com/Kynema/kynema-sgf/blob/development/unit_tests/wind_energy/test_abl_init.cpp>`_

   This is an example that uses :class:`MeshTest` to generate a test mesh and
   test the ABL initial conditions generator algorithms.

`test_refinement.cpp <https://github.com/Kynema/kynema-sgf/blob/development/unit_tests/utilities/test_refinement.cpp>`_

   This is an advanced example that test the user-defined nested mesh refinement
   algorithm by creating a test fixture that is capable of adaptive mesh
   refinement based on the criteria.
