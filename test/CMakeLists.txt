#=============================================================================
# Kynema-SGF Testing
#=============================================================================

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake)

if(KYNEMA_SGF_ENABLE_MASA AND NOT KYNEMA_SGF_ENABLE_MPI)
  message(WARNING "Running verification tests without MPI enabled will require long run times")
endif()

if(KYNEMA_SGF_TEST_WITH_FCOMPARE)
  if("${KYNEMA_SGF_REFERENCE_GOLDS_DIRECTORY}" STREQUAL "")
    message(FATAL_ERROR "To reference gold files, KYNEMA_SGF_REFERENCE_GOLDS_DIRECTORY must be set and exist")
  else()
    set(GOLD_FILES_DIRECTORY ${KYNEMA_SGF_REFERENCE_GOLDS_DIRECTORY}/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})
    message(STATUS "Test golds directory for fcompare: ${GOLD_FILES_DIRECTORY}")
  endif()
endif()

if(KYNEMA_SGF_SAVE_GOLDS)
  if("${KYNEMA_SGF_SAVED_GOLDS_DIRECTORY}" STREQUAL "")
    message(FATAL_ERROR "To save gold files, KYNEMA_SGF_SAVED_GOLDS_DIRECTORY must be set and the directory exist")
  else()
    if(EXISTS ${KYNEMA_SGF_SAVED_GOLDS_DIRECTORY})
      set(SAVED_GOLDS_DIRECTORY ${KYNEMA_SGF_SAVED_GOLDS_DIRECTORY}/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})
      message(STATUS "Gold files will be saved to: ${SAVED_GOLDS_DIRECTORY}")
    else()
      message(FATAL_ERROR "Specified directory for saving gold files does not exist: ${KYNEMA_SGF_SAVED_GOLDS_DIRECTORY}")
    endif()
  endif()
endif()

# Have CMake discover the number of cores on the node
include(ProcessorCount)
ProcessorCount(PROCESSES)

#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================
macro(setup_test)
    set(CURRENT_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${TEST_NAME})
    set(CURRENT_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${TEST_NAME})
    set(PLOT_GOLD ${GOLD_FILES_DIRECTORY}/${TEST_NAME}/plt00006)
    set(PLOT_TEST ${CURRENT_TEST_BINARY_DIR}/plt00006)
    if("${TEST_NAME}" MATCHES "hdf5$")
      set(PLOT_TEST ${PLOT_TEST}.h5)
    endif()
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    set(RUNTIME_OPTIONS "time.max_step=6 io.plot_file=plt time.plot_interval=6 amrex.the_arena_is_managed=0 amrex.abort_on_out_of_gpu_memory=1 amrex.the_arena_init_size=0")
    if((NOT KYNEMA_SGF_TEST_WITH_FCOMPARE) AND (NOT KYNEMA_SGF_SAVE_GOLDS)) #Avoid all I/O otherwise
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} time.plot_interval=-1 time.checkpoint_interval=-1")
    endif()
    if(KYNEMA_SGF_ENABLE_FPE_TRAP_FOR_TESTS)
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} amrex.signal_handling=1 amrex.fpe_trap_invalid=1 amrex.fpe_trap_zero=1 amrex.fpe_trap_overflow=1")
    else()
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} amrex.signal_handling=0")
    endif()
    if(KYNEMA_SGF_ENABLE_MPI)
      if(KYNEMA_SGF_ENABLE_CUDA)
        set(TEST_NP 2)
      else()
        set(TEST_NP 4)
      endif()
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${TEST_NP} ${MPIEXEC_PREFLAGS}")
    else()
      set(TEST_NP 1)
      unset(MPI_COMMANDS)
    endif()
    if((KYNEMA_SGF_ENABLE_CUDA OR KYNEMA_SGF_ENABLE_ROCM OR KYNEMA_SGF_ENABLE_SYCL) OR (KYNEMA_SGF_ENABLE_HYPRE AND (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")))
      set(FCOMPARE_TOLERANCE "-r 1e-10 --abs_tol 1.0e-12")
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} nodal_proj.deterministic=1 mac_proj.deterministic=1 diffusion.deterministic=1 io.skip_outputs=p gp")
    endif()
    if(KYNEMA_SGF_SAVE_GOLDS)
      file(MAKE_DIRECTORY ${SAVED_GOLDS_DIRECTORY}/${TEST_NAME})
      set(SAVE_GOLDS_COMMAND "&& cp -R ${PLOT_TEST} ${SAVED_GOLDS_DIRECTORY}/${TEST_NAME}/")
    endif()
    if(KYNEMA_SGF_TEST_WITH_FCOMPARE AND (NOT "${TEST_NAME}" MATCHES "hdf5$"))
      set(FCOMPARE_COMMAND "&& CUDA_LAUNCH_BLOCKING=1 ${MPI_COMMANDS} ${FCOMPARE_EXE} ${MPIEXEC_POSTFLAGS} ${FCOMPARE_TOLERANCE} ${PLOT_GOLD} ${PLOT_TEST}")
    endif()
endmacro(setup_test)

function(set_test_cost_if_provided TEST_NAME)
    if(ARGC GREATER 1)
      set_tests_properties(${TEST_NAME} PROPERTIES COST "${ARGV1}")
    endif()
endfunction(set_test_cost_if_provided)

# Standard regression test
function(add_test_r TEST_NAME)
    setup_test()
    add_test(${TEST_NAME} bash -c "set -o pipefail && ${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${kynema_sgf_exe_name} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.inp ${RUNTIME_OPTIONS} 2>&1 | tee ${TEST_NAME}.log ${SAVE_GOLDS_COMMAND} ${FCOMPARE_COMMAND}")
    # Set properties for test
    set_tests_properties(${TEST_NAME} PROPERTIES
                         TIMEOUT 14400
                         PROCESSORS ${TEST_NP}
                         WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
                         LABELS "regression"
                         ATTACHED_FILES "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
                         ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/Backtrace.0")
    set_test_cost_if_provided(${TEST_NAME} ${ARGN})
endfunction(add_test_r)

# Regression tests included in debug CI
function(add_test_rd TEST_NAME)
    add_test_r(${TEST_NAME} ${ARGN})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "regression;debug_ci")
endfunction(add_test_rd)

# Regression tests excluded from CI
function(add_test_re TEST_NAME)
    add_test_r(${TEST_NAME} ${ARGN})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "regression;no_ci")
endfunction(add_test_re)

# Regression test and excluded from CI with dependency
function(add_test_red TEST_NAME TEST_DEPENDENCY)
    if(ARGC GREATER 2)
      add_test_re(${TEST_NAME} ${ARGV2})
    else()
      add_test_re(${TEST_NAME})
    endif()
    set_tests_properties(${TEST_NAME} PROPERTIES FIXTURES_REQUIRED fixture_${TEST_DEPENDENCY})
    set_tests_properties(${TEST_DEPENDENCY} PROPERTIES FIXTURES_SETUP fixture_${TEST_DEPENDENCY})
endfunction(add_test_red)

# Run the input of SOURCE_TEST with KLAxell and with KLAxellSeparation and
# compare the two plot files with each other (no gold files). Further
# arguments are runtime options for the KLAxellSeparation run; options after
# COMMON_OPTIONS are given to both runs. By default the two plot files must
# agree; with RELATIVE_TOLERANCE <value> and ABSOLUTE_TOLERANCE <value> each
# variable must agree to either tolerance (for numerical options that keep the
# steady state; the absolute tolerance covers variables that are zero up to
# round-off), and with the EXPECT_DIFFERENT argument the KLAxellSeparation plot
# file must be free of NaNs and differ from KLAxell.
function(add_test_klaxell_pair TEST_NAME SOURCE_TEST)
    set(PAIR_NAME ${TEST_NAME})
    set(TEST_NAME ${SOURCE_TEST})
    setup_test()
    set(PAIR_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/${PAIR_NAME})
    cmake_parse_arguments(PAIR "EXPECT_DIFFERENT" "RELATIVE_TOLERANCE;ABSOLUTE_TOLERANCE" "COMMON_OPTIONS" ${ARGN})
    set(EXPECT_DIFFERENT ${PAIR_EXPECT_DIFFERENT})
    list(JOIN PAIR_UNPARSED_ARGUMENTS " " SEPARATION_OPTIONS)
    list(JOIN PAIR_COMMON_OPTIONS " " COMMON_OPTIONS)
    if(DEFINED PAIR_RELATIVE_TOLERANCE)
      set(FCOMPARE_TOLERANCE "${FCOMPARE_TOLERANCE} -r ${PAIR_RELATIVE_TOLERANCE}")
    endif()
    if(DEFINED PAIR_ABSOLUTE_TOLERANCE)
      set(FCOMPARE_TOLERANCE "${FCOMPARE_TOLERANCE} --abs_tol ${PAIR_ABSOLUTE_TOLERANCE}")
    endif()
    unset(PAIR_COMMAND)
    foreach(MODEL IN ITEMS KLAxell KLAxellSeparation)
      file(MAKE_DIRECTORY ${PAIR_DIR}/${MODEL})
      file(COPY ${TEST_FILES} DESTINATION "${PAIR_DIR}/${MODEL}/")
      set(MODEL_OPTIONS "turbulence.model=${MODEL} ${COMMON_OPTIONS}")
      if("${MODEL}" STREQUAL "KLAxellSeparation")
        string(APPEND MODEL_OPTIONS " ${SEPARATION_OPTIONS}")
      endif()
      string(APPEND PAIR_COMMAND "cd ${PAIR_DIR}/${MODEL} && ${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${kynema_sgf_exe_name} ${MPIEXEC_POSTFLAGS} ${PAIR_DIR}/${MODEL}/${SOURCE_TEST}.inp ${RUNTIME_OPTIONS} ${MODEL_OPTIONS} 2>&1 | tee ${PAIR_NAME}_${MODEL}.log && ")
    endforeach()
    set(COMPARE_COMMAND "CUDA_LAUNCH_BLOCKING=1 ${MPI_COMMANDS} ${FCOMPARE_EXE} ${MPIEXEC_POSTFLAGS} ${FCOMPARE_TOLERANCE} ${PAIR_DIR}/KLAxell/plt00006 ${PAIR_DIR}/KLAxellSeparation/plt00006")
    if(EXPECT_DIFFERENT)
      string(REPLACE "amrex_fcompare" "amrex_fnan" FNAN_EXE "${FCOMPARE_EXE}")
      string(APPEND PAIR_COMMAND "${FNAN_EXE} ${PAIR_DIR}/KLAxellSeparation/plt00006 && ! ${COMPARE_COMMAND}")
    else()
      string(APPEND PAIR_COMMAND "${COMPARE_COMMAND}")
    endif()
    add_test(${PAIR_NAME} bash -c "set -o pipefail && ${PAIR_COMMAND}")
    set_tests_properties(${PAIR_NAME} PROPERTIES
                         TIMEOUT 14400
                         PROCESSORS ${TEST_NP}
                         WORKING_DIRECTORY "${PAIR_DIR}/"
                         LABELS "regression"
                         ATTACHED_FILES "${PAIR_DIR}/KLAxellSeparation/${PAIR_NAME}_KLAxellSeparation.log")
endfunction(add_test_klaxell_pair)

# Verification test using multiple resolutions
function(add_test_v TEST_NAME LIST_OF_GRID_SIZES)
    setup_test()
    unset(MASTER_RUN_COMMAND)
    # Get last item in resolution list so we can find out when we are on the last item in our loop
    list(GET LIST_OF_GRID_SIZES -1 LAST_GRID_SIZE_IN_LIST)
    foreach(GRID_SIZE IN LISTS LIST_OF_GRID_SIZES)
      file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE})
      file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
      file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}/")
      set(NCELLS "${GRID_SIZE} ${GRID_SIZE} ${GRID_SIZE}")
      set(RUN_COMMAND_${GRID_SIZE} "${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${kynema_sgf_exe_name} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}/${TEST_NAME}.inp")
      set(RUNTIME_OPTIONS_${GRID_SIZE} "amrex.throw_exception=1 amrex.signal_handling=0 amr.n_cell=${NCELLS}")
      string(APPEND MASTER_RUN_COMMAND "cd ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}")
      string(APPEND MASTER_RUN_COMMAND " && ")
      string(APPEND MASTER_RUN_COMMAND "${RUN_COMMAND_${GRID_SIZE}} ${RUNTIME_OPTIONS_${GRID_SIZE}} 2>&1 | tee ${TEST_NAME}_${GRID_SIZE}.log")
      # Add another " && " unless we are on the last resolution in the list
      if(NOT ${GRID_SIZE} EQUAL ${LAST_GRID_SIZE_IN_LIST})
        string(APPEND MASTER_RUN_COMMAND " && ")
      endif()
    endforeach()
    list(JOIN LIST_OF_GRID_SIZES " " STRING_OF_GRID_SIZES)
    add_test(${TEST_NAME} bash -c "set -o pipefail && ${MASTER_RUN_COMMAND} && cd ${CURRENT_TEST_BINARY_DIR} && ${PYTHON_EXECUTABLE} ${CURRENT_TEST_SOURCE_DIR}/plotter.py -f ${STRING_OF_GRID_SIZES}")
    set_tests_properties(${TEST_NAME} PROPERTIES TIMEOUT 14400 PROCESSORS ${TEST_NP} WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}" LABELS "verification;no_ci" ATTACHED_FILES "${CURRENT_TEST_BINARY_DIR}/plots.pdf")
    set_test_cost_if_provided(${TEST_NAME} ${ARGN})
endfunction(add_test_v)

# Standard unit test
function(add_test_u TEST_NAME)
    setup_test()
    set(TEST_NP 1)
    if(KYNEMA_SGF_ENABLE_MPI)
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${TEST_NP} ${MPIEXEC_PREFLAGS}")
    else()
      unset(MPI_COMMANDS)
    endif()
    add_test(${TEST_NAME} bash -c "${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${kynema_sgf_unit_test_exe_name}")
    set_tests_properties(${TEST_NAME} PROPERTIES
                         TIMEOUT 500
                         PROCESSORS ${TEST_NP}
                         WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
                         LABELS "unit")
    set_test_cost_if_provided(${TEST_NAME} ${ARGN})
endfunction(add_test_u)

# Add test that also checks post-processing scripts: sampling
function(add_test_pps TEST_NAME FORMAT NINT NSTEPS)
    setup_test()
    set(ADDL_RUNTIME_OPTIONS "sampling.output_format=${FORMAT} sampling.output_interval=${NINT} time.plot_interval=-1 time.checkpoint_interval=-1 time.max_step=${NSTEPS}")
    add_test("${TEST_NAME}_sampling_${FORMAT}" bash -c "set -o pipefail && ${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${kynema_sgf_exe_name} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.inp ${RUNTIME_OPTIONS} ${ADDL_RUNTIME_OPTIONS} 2>&1 | tee ${TEST_NAME}_sampling_${FORMAT}.log && ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/sampling_${TEST_NAME}_${FORMAT}.py")
    # Set properties for test
    set_tests_properties("${TEST_NAME}_sampling_${FORMAT}" PROPERTIES
                          TIMEOUT 14400
                          PROCESSORS ${TEST_NP}
                          WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
                          LABELS "post_processing;no_ci"
                          ATTACHED_FILES "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}_sampling_${FORMAT}.log;plot_sampling_${FORMAT}.pdf")
    set_test_cost_if_provided("${TEST_NAME}_sampling_${FORMAT}" ${ARGN})
endfunction(add_test_pps)

#=============================================================================
# Unit tests
#=============================================================================
add_test_u(unit_tests 122.19)

#=============================================================================
# Subset of regression tests for both debug and release CI
#=============================================================================
add_test_rd(abl_godunov 483.65)
add_test_rd(abl_godunov_timetable 488.01)
add_test_rd(abl_godunov_static_refinement 568.44)
add_test_rd(abl_multiphase_laminar 2544.41)
add_test_rd(abl_sampling 727.79)
add_test_rd(abl_specifiedmol_neutral 565.48)
add_test_rd(act_fixed_wing 334.44)
add_test_rd(boussinesq_bubble_godunov 307.55)
add_test_rd(channel_builder_single_phase_periodic 357.76)
add_test_rd(channel_godunov_laminar 19.24)
add_test_rd(channel_kwsstiddes 200.71)
add_test_rd(channel_smagorinsky_analytical 11.85)
add_test_rd(channel_terrain_walls 17.91)
add_test_rd(dam_break_godunov 837.87)
add_test_rd(ekman_spiral 19.04)
add_test_rd(forest_drag 570.78)
add_test_rd(forest_drag_point_cloud 581.64)
add_test_rd(freestream_godunov_inout_refine2 36.6)
add_test_rd(halfchannel_symmetricwall 8.77)
add_test_rd(halfchannel_zerogradient 8.1)
add_test_rd(hbl_godunov 706.03)
add_test_rd(joukowsky_disk 805.89)
add_test_rd(ow_current_bathymetry 401.47)
add_test_rd(ow_current_wing 1770.65)
add_test_rd(rain_drop 1278.15)
add_test_rd(terrain_box 1065.55)
add_test_rd(udf_refinement 1125.3)
add_test_rd(uniform_ct_disk_dynamic_adaptation 891.22)
add_test_rd(vortex_patch_godunov 154.76)


#=============================================================================
# Regression tests (release CI)
#=============================================================================
add_test_r(abl_amd_wenoz 639.03)
add_test_r(abl_anelastic 190.32)
add_test_r(abl_bndry_output_amr_native 331.45)
add_test_r(abl_bndry_output_native 457.36)
add_test_r(abl_godunov_geostrophic_timetable 505.79)
add_test_r(abl_godunov_mpl 894.46)
add_test_r(abl_godunov_mpl_amr 394.34)
add_test_r(abl_godunov_rayleigh_damping 2166.2)
add_test_r(abl_godunov_scalar_velocity_solve 450.92)
add_test_r(abl_godunov_segregated_velocity_solve 444.48)
add_test_r(abl_godunov_static_refinement_rr4 1252.03)
add_test_r(abl_godunov_stats_max_level 503.99)
add_test_r(abl_godunov_wenoz_fixedpt 931.4)
add_test_r(abl_ksgsm84_godunov 465.52)
add_test_r(abl_rans_convergence)
add_test_r(abl_multiphase_laminar_2planes 2488.54)
add_test_r(abl_specifiedmol_stable 487.57)
add_test_r(abl_specifiedmol_unstable 488.46)
if(KYNEMA_SGF_TEST_WITH_FCOMPARE)
  add_test_klaxell_pair(klaxellseparation_neutral abl_specifiedmol_neutral)
  add_test_klaxell_pair(klaxellseparation_stable abl_specifiedmol_stable)
  add_test_klaxell_pair(klaxellseparation_unstable abl_specifiedmol_unstable)
  add_test_klaxell_pair(klaxellseparation_hill abl_klaxell_hill)
  add_test_klaxell_pair(klaxellseparation_stable_sensor abl_specifiedmol_stable KLAxellSeparation.pressure_gradient_sensor=true)
  add_test_klaxell_pair(klaxellseparation_hill_sensor abl_klaxell_hill KLAxellSeparation.pressure_gradient_sensor=true)
  add_test_klaxell_pair(klaxellseparation_neutral_realizable_cmu abl_specifiedmol_neutral KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.realizable_cmu=true KLAxellSeparation_coeffs.gate_relaxation_time=0)
  add_test_klaxell_pair(klaxellseparation_stable_realizable_cmu abl_specifiedmol_stable KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.realizable_cmu=true KLAxellSeparation_coeffs.gate_relaxation_time=0)
  add_test_klaxell_pair(klaxellseparation_hill_realizable_cmu_off abl_klaxell_hill KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.realizable_cmu=true KLAxellSeparation_coeffs.realizable_cmu_strength=0)
  add_test_klaxell_pair(klaxellseparation_hill_realizable_cmu abl_klaxell_hill EXPECT_DIFFERENT KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.realizable_cmu=true KLAxellSeparation_coeffs.gate_relaxation_time=0)
  add_test_klaxell_pair(klaxellseparation_neutral_velocity_sensor abl_specifiedmol_neutral KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.sensor_source=velocity KLAxellSeparation.realizable_cmu=true)
  add_test_klaxell_pair(klaxellseparation_stable_velocity_sensor abl_specifiedmol_stable KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.sensor_source=velocity KLAxellSeparation.realizable_cmu=true)
  add_test_klaxell_pair(klaxellseparation_hill_velocity_sensor abl_klaxell_hill KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.sensor_source=velocity)
  add_test_klaxell_pair(klaxellseparation_hill_velocity_realizable_cmu abl_klaxell_hill EXPECT_DIFFERENT KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.sensor_source=velocity KLAxellSeparation.realizable_cmu=true KLAxellSeparation_coeffs.gate_relaxation_time=0)
  add_test_klaxell_pair(klaxellseparation_neutral_gate_relaxation abl_specifiedmol_neutral KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.realizable_cmu=true KLAxellSeparation_coeffs.gate_relaxation_time=10.0)
  add_test_klaxell_pair(klaxellseparation_stable_gate_relaxation abl_specifiedmol_stable KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.realizable_cmu=true KLAxellSeparation_coeffs.gate_relaxation_time=10.0)
  add_test_klaxell_pair(klaxellseparation_hill_gate_relaxation_diagnostic abl_klaxell_hill KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation_coeffs.gate_relaxation_time=10.0)
  add_test_klaxell_pair(klaxellseparation_hill_gate_relaxation abl_klaxell_hill EXPECT_DIFFERENT KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.realizable_cmu=true KLAxellSeparation_coeffs.gate_relaxation_time=10.0)
  add_test_klaxell_pair(klaxellseparation_neutral_production_cap abl_specifiedmol_neutral KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.production_cap=true KLAxellSeparation_coeffs.gate_relaxation_time=0)
  add_test_klaxell_pair(klaxellseparation_stable_production_cap abl_specifiedmol_stable KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.production_cap=true KLAxellSeparation_coeffs.gate_relaxation_time=0)
  add_test_klaxell_pair(klaxellseparation_hill_production_cap_off abl_klaxell_hill KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.production_cap=true KLAxellSeparation_coeffs.sensor_threshold=1000)
  add_test_klaxell_pair(klaxellseparation_hill_production_cap abl_klaxell_hill EXPECT_DIFFERENT KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.production_cap=true KLAxellSeparation_coeffs.production_cap_ratio=1 KLAxellSeparation_coeffs.gate_relaxation_time=0)
  add_test_klaxell_pair(klaxellseparation_neutral_production_cap_gate_relaxation abl_specifiedmol_neutral KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.realizable_cmu=true KLAxellSeparation.production_cap=true KLAxellSeparation_coeffs.gate_relaxation_time=10.0)
  add_test_klaxell_pair(klaxellseparation_hill_production_cap_gate_relaxation abl_klaxell_hill EXPECT_DIFFERENT KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.production_cap=true KLAxellSeparation_coeffs.production_cap_ratio=1 KLAxellSeparation_coeffs.gate_relaxation_time=10.0)
  add_test_klaxell_pair(klaxellseparation_neutral_destruction_boost abl_specifiedmol_neutral KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.destruction_boost=true KLAxellSeparation_coeffs.gate_relaxation_time=0)
  add_test_klaxell_pair(klaxellseparation_stable_destruction_boost abl_specifiedmol_stable KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.destruction_boost=true KLAxellSeparation_coeffs.gate_relaxation_time=0)
  add_test_klaxell_pair(klaxellseparation_hill_destruction_boost_off abl_klaxell_hill KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.destruction_boost=true KLAxellSeparation_coeffs.sensor_threshold=1000)
  add_test_klaxell_pair(klaxellseparation_hill_destruction_boost abl_klaxell_hill EXPECT_DIFFERENT KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.destruction_boost=true KLAxellSeparation_coeffs.destruction_boost_factor=2 KLAxellSeparation_coeffs.gate_relaxation_time=0)
  add_test_klaxell_pair(klaxellseparation_hill_destruction_boost_gate_relaxation abl_klaxell_hill EXPECT_DIFFERENT KLAxellSeparation.pressure_gradient_sensor=true KLAxellSeparation.destruction_boost=true KLAxellSeparation_coeffs.destruction_boost_factor=2 KLAxellSeparation_coeffs.gate_relaxation_time=10.0)
  add_test_klaxell_pair(klaxellseparation_neutral_curvature_richardson abl_specifiedmol_neutral KLAxellSeparation.curvature_correction=true KLAxellSeparation.curvature_model=richardson)
  add_test_klaxell_pair(klaxellseparation_stable_curvature_richardson abl_specifiedmol_stable KLAxellSeparation.curvature_correction=true KLAxellSeparation.curvature_model=richardson)
  add_test_klaxell_pair(klaxellseparation_hill_curvature_richardson abl_klaxell_hill EXPECT_DIFFERENT KLAxellSeparation.curvature_correction=true KLAxellSeparation.curvature_model=richardson)
  add_test_klaxell_pair(klaxellseparation_neutral_curvature_rotation abl_specifiedmol_neutral KLAxellSeparation.curvature_correction=true KLAxellSeparation.curvature_model=rotation_function)
  add_test_klaxell_pair(klaxellseparation_stable_curvature_rotation abl_specifiedmol_stable KLAxellSeparation.curvature_correction=true KLAxellSeparation.curvature_model=rotation_function)
  add_test_klaxell_pair(klaxellseparation_hill_curvature_rotation abl_klaxell_hill EXPECT_DIFFERENT KLAxellSeparation.curvature_correction=true KLAxellSeparation.curvature_model=rotation_function)
  add_test_klaxell_pair(klaxellseparation_neutral_implicit_dissipation abl_specifiedmol_neutral RELATIVE_TOLERANCE 5.0e-3 ABSOLUTE_TOLERANCE 1.0e-10 KLAxellSeparation.implicit_dissipation=true)
  add_test_klaxell_pair(klaxellseparation_stable_implicit_dissipation abl_specifiedmol_stable RELATIVE_TOLERANCE 5.0e-3 ABSOLUTE_TOLERANCE 1.0e-10 KLAxellSeparation.implicit_dissipation=true)
  add_test_klaxell_pair(klaxellseparation_neutral_implicit_dissipation_crank_nicolson abl_specifiedmol_neutral RELATIVE_TOLERANCE 5.0e-3 ABSOLUTE_TOLERANCE 1.0e-10 COMMON_OPTIONS incflo.diffusion_type=1 KLAxellSeparation.implicit_dissipation=true)
  add_test_klaxell_pair(klaxellseparation_hill_implicit_dissipation abl_klaxell_hill EXPECT_DIFFERENT KLAxellSeparation.implicit_dissipation=true)
endif()
add_test_r(abl_stable 426.56)
add_test_r(abl_stable_different_roughness 424.42)
add_test_r(abl_subvolume 623.01)
add_test_r(abl_surf_temp_timetable 200.9)
add_test_r(abl_unstable 161.84)
add_test_r(abl_unstable_constant_wall_model 155.02)
add_test_r(abl_unstable_local_wall_model 155.56)
add_test_r(abl_unstable_schumann_wall_model 149.02)
add_test_r(abl_waves_terrain 945.4)
add_test_r(act_abl_joukowskydisk 1402.12)
add_test_r(act_abl_uniformctdisk 1334.13)
add_test_r(act_fixed_wing_fllc 267.08)
add_test_r(act_flat_plate 319.3)
add_test_r(act_moving_wing 111.98)
add_test_r(act_pitching_wing_2D 425.74)
add_test_r(act_sector 219.22)
add_test_r(act_sector_quad 488.73)
add_test_r(act_drone_quad 496.13)
add_test_r(box_refinement 127.92)
add_test_r(burggraf_flow 135.37)
add_test_r(channel_builder_multiphase_contraction 687.69)
add_test_r(channel_kwsst 170.16)
add_test_r(channel_kwsst_sust 170.8)
add_test_r(ctv_godunov_plm 95.44)
add_test_r(cylinder_refinement 1476.6)
add_test_r(fat_cored_vortex_ring 252.68)
add_test_r(freestream_godunov 59.67)
add_test_r(freestream_godunov_inout 14.41)
add_test_r(freestream_godunov_inout_refine1 33.65)
add_test_r(inertial_drop 1162.25)
add_test_r(inflow_godunov_amr 148.06)
add_test_r(joukowsky_disk_awc 709.76)
add_test_r(linear_godunov_amr 456.78)
add_test_r(scalar_advection_refined 51.86)
add_test_r(scalar_advection_uniform 10.13)
add_test_r(sloshing_tank 974.62)
add_test_r(terrain_box_amr 2550.02)
add_test_r(tgv_godunov 586.68)
add_test_r(tgv_godunov_plm 592.59)
add_test_r(uniform_ct_disk 649.68)
add_test_r(uniform_ct_disk_gaussian 651.13)
add_test_r(vortex_dipole_inout_refine 47.87)
add_test_r(vortex_dipole_wall_collision 729.89)
add_test_r(vortex_patch_scalar_vel 26.63)
add_test_r(vortex_ring_collision 90.14)
add_test_r(zalesak_disk_godunov 81.8)
add_test_r(zalesak_disk_scalar_vel 15.32)


#=============================================================================
# Regression tests excluded from CI
#=============================================================================
add_test_re(abl_bds 1454.23)
add_test_re(abl_godunov_cn 402.17)
add_test_re(abl_godunov_explicit 317.68)
add_test_re(abl_godunov_nolim 393.33)
add_test_re(abl_godunov_noncons 400.43)
add_test_re(abl_godunov_plm 383.42)
add_test_re(abl_godunov_ppm 366.99)
add_test_re(abl_godunov_weno 369.86)
add_test_re(abl_mol 451.44)
add_test_re(abl_mol_cn 468.63)
add_test_re(abl_mol_explicit 328.06)
add_test_re(boussinesq_bubble_mol 263.12)
add_test_re(ctv_bds 271.6)
add_test_re(ctv_mol_mesh_map 479.77)
add_test_re(ctv_mol_mesh_map_explicit 425.39)
add_test_re(freestream_bds 175.26)
add_test_re(ib_ctv_godunov_weno 82.27)
add_test_re(ib_cylinder_Re_300 994.38)
add_test_re(ib_sphere_Re_100 1253.34)
add_test_re(inflow_bds_amr 906.66)
add_test_re(linear_bds_amr 965.48)
add_test_re(rayleigh_taylor_godunov 808.09)
add_test_re(rayleigh_taylor_mol 1186.89)
add_test_re(tgv_mol 634.63)

if("${KYNEMA_SGF_PRECISION}" STREQUAL "DOUBLE")
  add_test_re(ow_linear 1282.28)
  add_test_re(ow_linear_init_waves 338.87)
  add_test_re(ow_stokes 1097.54)
  add_test_re(rankine 101.32)
  add_test_re(rankine-sym 98.3)
endif()

if(KYNEMA_SGF_ENABLE_NETCDF)
  add_test_re(abl_bndry_output 476.96)
  add_test_re(abl_bndry_output_amr_inflow 341)
  add_test_re(abl_bndry_output_amr_upper 285.1)
  add_test_re(abl_kosovic_neutral 485.44)
  add_test_re(abl_kosovic_neutral_ib 1860.11)
  add_test_re(abl_meso_input_dpa 62.36)
  add_test_re(abl_meso_input_ipa 58.95)
  add_test_re(abl_meso_tendency 54.48)
  add_test_re(abl_sampling_netcdf 557.72)
  add_test_re(abl_virtual_lidar 370.84)
  add_test_re(abl_wallrans_neutral 781.92)
  add_test_re(nrel_precursor 715.27)
endif()

if(KYNEMA_SGF_ENABLE_MASA)
  add_test_re(mms_bds)
  add_test_re(mms_godunov)
  add_test_re(mms_godunov_plm)
  add_test_re(mms_mol)
endif()

if(KYNEMA_SGF_ENABLE_HYPRE)
  add_test_re(abl_godunov_hypre 350.52)
  add_test_re(channel_kwsst_hypre 641.67)
  add_test_re(channel_mol_mesh_map_x 37.42)
  add_test_re(channel_mol_mesh_map_x_seg_vel_solve 24.21)
  add_test_re(channel_mol_mesh_map_y 34.98)
  add_test_re(channel_mol_mesh_map_z 36.08)
endif()

if(KYNEMA_SGF_ENABLE_HDF5)
  add_test_re(abl_stable_hdf5 372.76)
  if(KYNEMA_SGF_ENABLE_HDF5_ZFP)
    add_test_re(abl_stable_zfp_hdf5 357.71)
  endif()
endif()

if(KYNEMA_SGF_ENABLE_ASCENT)
  add_test_re(abl_godunov_ascent)
endif()

if(KYNEMA_SGF_ENABLE_W2A)
  add_test_re(abl_multiphase_w2a 1486.31)
  add_test_re(ow_w2a 1612.5)
  add_test_re(ow_w2a_nonperiodic 243.03)
  add_test_re(ow_w2a_nwt_2d 478.31)
  add_test_re(ow_w2a_nwt_3d 455.64)
  if(KYNEMA_SGF_TEST_WITH_FCOMPARE)
    add_test_red(abl_w2a_terrain abl_bndry_output_native 471.55)
  endif()
endif()

if(KYNEMA_SGF_ENABLE_OPENFAST AND KYNEMA_SGF_OPENFAST_VERSION VERSION_LESS "4")
  set(ACT_UNIFORM_ALM_TEST_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/test_files/act_uniform_alm)
  set(ACT_UNIFORM_ALM_TEST_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/test_files/act_uniform_alm)
  set(RTEST_ACT_UNIFORM_ALM_TEST_BINARY_DIR ${ACT_UNIFORM_ALM_TEST_BINARY_DIR}/r-test)
  set(RTEST_TAG "v3.5.4")
  set(5MW_BASELINE_DIR glue-codes/openfast/5MW_Baseline)
  set(5MW_LAND_DIR glue-codes/openfast-cpp/5MW_Land_DLL_WTurb_cpp)
  set(AERO_FILE NRELOffshrBsline5MW_Onshore_AeroDyn15.dat)
  set(ELASTO_FILE NRELOffshrBsline5MW_Onshore_ElastoDyn.dat)
  set(SERVO_FILE NRELOffshrBsline5MW_Onshore_ServoDyn.dat)
  file(MAKE_DIRECTORY ${ACT_UNIFORM_ALM_TEST_BINARY_DIR})
  if(NOT EXISTS "${RTEST_ACT_UNIFORM_ALM_TEST_BINARY_DIR}")
    execute_process(
      COMMAND bash -c "${GIT_EXECUTABLE} clone -n --depth=1 --filter=tree:0 https://github.com/OpenFAST/r-test.git"
      WORKING_DIRECTORY
      "${ACT_UNIFORM_ALM_TEST_BINARY_DIR}"
      COMMAND_ERROR_IS_FATAL ANY)
    execute_process(
      COMMAND bash -c "${GIT_EXECUTABLE} fetch --all --tags"
      WORKING_DIRECTORY
      "${RTEST_ACT_UNIFORM_ALM_TEST_BINARY_DIR}"
      COMMAND_ERROR_IS_FATAL ANY)
    execute_process(
      COMMAND bash -c "${GIT_EXECUTABLE} sparse-checkout set --no-cone ${5MW_BASELINE_DIR} ${5MW_LAND_DIR}"
      WORKING_DIRECTORY
      "${RTEST_ACT_UNIFORM_ALM_TEST_BINARY_DIR}"
      COMMAND_ERROR_IS_FATAL ANY)
    execute_process(
      COMMAND bash -c "${GIT_EXECUTABLE} checkout ${RTEST_TAG}"
      WORKING_DIRECTORY
      "${RTEST_ACT_UNIFORM_ALM_TEST_BINARY_DIR}"
      COMMAND_ERROR_IS_FATAL ANY)
    execute_process(
      COMMAND bash -c "${GIT_EXECUTABLE} apply ${ACT_UNIFORM_ALM_TEST_SOURCE_DIR}/custom_openfast_options.patch"
      WORKING_DIRECTORY
      "${RTEST_ACT_UNIFORM_ALM_TEST_BINARY_DIR}"
      COMMAND_ERROR_IS_FATAL ANY)
  endif()
  add_test_re(act_uniform_alm)
  if(KYNEMA_SGF_ENABLE_NETCDF AND KYNEMA_SGF_TEST_WITH_FCOMPARE)
    add_test_red(act_uniform_alm_restart act_uniform_alm)
  endif()
endif()

if(KYNEMA_SGF_ENABLE_KYNEMA_FMB)
  add_test_re(act_kynema_fmb_alm)
  add_test_re(act_kynema_fmb_adm)
  add_test_red(act_kynema_fmb_alm_restart act_kynema_fmb_alm)
endif()

if(KYNEMA_SGF_ENABLE_FFT)
  add_test_re(abl_godunov_nofft)
endif()
#=============================================================================
# Regression tests excluded from CI with a test dependency
#=============================================================================
if(KYNEMA_SGF_TEST_WITH_FCOMPARE)
  add_test_red(abl_bndry_input_native abl_bndry_output_native 116.58)
  add_test_red(abl_bndry_input_native_inout abl_bndry_output_native 182.49)
  add_test_red(abl_godunov_restart abl_godunov 147.85)
  add_test_red(abl_bndry_input_amr_native abl_bndry_output_native 317.8)
  add_test_red(abl_bndry_input_amr_native_xhi abl_bndry_output_native 315.98)
  add_test_red(abl_bndry_input_amr_native_mlbc abl_bndry_output_amr_native 74.55)
  add_test_red(abl_godunov_forcetimetable abl_godunov_timetable 131.43)
  add_test_red(abl_multiphase_laminar_input abl_multiphase_laminar 635.14)
  add_test_red(abl_multiphase_laminar_inout abl_multiphase_laminar_2planes 638.39)
  add_test_red(channel_builder_jet_input channel_builder_multiphase_contraction 491.55)
  if(KYNEMA_SGF_ENABLE_NETCDF)
    add_test_red(abl_bndry_input abl_bndry_output 114.44)
    add_test_red(abl_bndry_input_init abl_bndry_output 328.43)
    add_test_red(abl_bndry_input_amr abl_bndry_output 330.97)
    add_test_red(abl_bndry_input_amr_inflow abl_bndry_output_amr_inflow 72.97)
    add_test_red(abl_bndry_input_amr_upper abl_bndry_output_amr_upper 62.08)
    add_test_red(nrel_terrain nrel_precursor 336.27)
    add_test_red(nrel_terrain_amr nrel_precursor 377.86)
  endif()
endif()

#=============================================================================
# Verification tests
#=============================================================================
if(KYNEMA_SGF_ENABLE_MASA)
  set(LIST_OF_GRID_SIZES 8 16 32 64)
  add_test_v(mms "${LIST_OF_GRID_SIZES}")
endif()

#=============================================================================
# Postprocessing tests
#=============================================================================
if (KYNEMA_SGF_TEST_WITH_PYTHON)
  add_test_pps(dam_break_godunov native 30 150)
  add_test_pps(dam_break_godunov ascii 30 150)
  if(KYNEMA_SGF_ENABLE_NETCDF)
    add_test_pps(dam_break_godunov netcdf 30 150)
  endif()
endif()
