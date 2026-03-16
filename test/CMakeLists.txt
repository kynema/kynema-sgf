#=============================================================================
# AMR-Wind Testing
#=============================================================================

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CTestCustom.cmake ${CMAKE_BINARY_DIR}/CTestCustom.cmake)

if(AMR_WIND_ENABLE_MASA AND NOT AMR_WIND_ENABLE_MPI)
  message(WARNING "Running verification tests without MPI enabled will require long run times")
endif()

if(AMR_WIND_TEST_WITH_FCOMPARE)
  if("${AMR_WIND_REFERENCE_GOLDS_DIRECTORY}" STREQUAL "")
    message(FATAL_ERROR "To reference gold files, AMR_WIND_REFERENCE_GOLDS_DIRECTORY must be set and exist")
  else()
    set(GOLD_FILES_DIRECTORY ${AMR_WIND_REFERENCE_GOLDS_DIRECTORY}/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})
    message(STATUS "Test golds directory for fcompare: ${GOLD_FILES_DIRECTORY}")
  endif()
endif()

if(AMR_WIND_SAVE_GOLDS)
  if("${AMR_WIND_SAVED_GOLDS_DIRECTORY}" STREQUAL "")
    message(FATAL_ERROR "To save gold files, AMR_WIND_SAVED_GOLDS_DIRECTORY must be set and the directory exist")
  else()
    if(EXISTS ${AMR_WIND_SAVED_GOLDS_DIRECTORY})
      set(SAVED_GOLDS_DIRECTORY ${AMR_WIND_SAVED_GOLDS_DIRECTORY}/${CMAKE_SYSTEM_NAME}/${CMAKE_CXX_COMPILER_ID}/${CMAKE_CXX_COMPILER_VERSION})
      message(STATUS "Gold files will be saved to: ${SAVED_GOLDS_DIRECTORY}")
    else()
      message(FATAL_ERROR "Specified directory for saving gold files does not exist: ${AMR_WIND_SAVED_GOLDS_DIRECTORY}")
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
    set(PLOT_GOLD ${GOLD_FILES_DIRECTORY}/${TEST_NAME}/plt00010)
    set(PLOT_TEST ${CURRENT_TEST_BINARY_DIR}/plt00010)
    if("${TEST_NAME}" MATCHES "hdf5$")
      set(PLOT_TEST ${PLOT_TEST}.h5)
    endif()
    file(MAKE_DIRECTORY ${CURRENT_TEST_BINARY_DIR})
    file(GLOB TEST_FILES "${CURRENT_TEST_SOURCE_DIR}/*")
    file(COPY ${TEST_FILES} DESTINATION "${CURRENT_TEST_BINARY_DIR}/")
    set(RUNTIME_OPTIONS "time.max_step=10 io.plot_file=plt time.plot_interval=10 amrex.the_arena_is_managed=0 amrex.abort_on_out_of_gpu_memory=1 amrex.the_arena_init_size=0")
    if((NOT AMR_WIND_TEST_WITH_FCOMPARE) AND (NOT AMR_WIND_SAVE_GOLDS)) #Avoid all I/O otherwise
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} time.plot_interval=-1 time.checkpoint_interval=-1")
    endif()
    if(AMR_WIND_ENABLE_FPE_TRAP_FOR_TESTS)
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} amrex.signal_handling=1 amrex.fpe_trap_invalid=1 amrex.fpe_trap_zero=1 amrex.fpe_trap_overflow=1")
    else()
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} amrex.signal_handling=0")
    endif()
    if(AMR_WIND_ENABLE_MPI)
      if(AMR_WIND_ENABLE_CUDA)
        set(TEST_NP 2)
      else()
        set(TEST_NP 4)
      endif()
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${TEST_NP} ${MPIEXEC_PREFLAGS}")
    else()
      set(TEST_NP 1)
      unset(MPI_COMMANDS)
    endif()
    if((AMR_WIND_ENABLE_CUDA OR AMR_WIND_ENABLE_ROCM OR AMR_WIND_ENABLE_SYCL) OR (AMR_WIND_ENABLE_HYPRE AND (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")))
      set(FCOMPARE_TOLERANCE "-r 1e-10 --abs_tol 1.0e-12")
      set(RUNTIME_OPTIONS "${RUNTIME_OPTIONS} io.skip_outputs=p gp")
    endif()
    if(AMR_WIND_SAVE_GOLDS)
      file(MAKE_DIRECTORY ${SAVED_GOLDS_DIRECTORY}/${TEST_NAME})
      set(SAVE_GOLDS_COMMAND "&& cp -R ${PLOT_TEST} ${SAVED_GOLDS_DIRECTORY}/${TEST_NAME}/")
    endif()
    if(AMR_WIND_TEST_WITH_FCOMPARE AND (NOT "${TEST_NAME}" MATCHES "hdf5$"))
      set(FCOMPARE_COMMAND "&& CUDA_LAUNCH_BLOCKING=1 ${MPI_COMMANDS} ${FCOMPARE_EXE} ${MPIEXEC_POSTFLAGS} ${FCOMPARE_TOLERANCE} ${PLOT_GOLD} ${PLOT_TEST}")
    endif()
endmacro(setup_test)

# Standard regression test
function(add_test_r TEST_NAME)
    setup_test()
    add_test(${TEST_NAME} bash -c "set -o pipefail && ${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${amr_wind_exe_name} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.inp ${RUNTIME_OPTIONS} 2>&1 | tee ${TEST_NAME}.log ${SAVE_GOLDS_COMMAND} ${FCOMPARE_COMMAND}")
    # Set properties for test
    set_tests_properties(${TEST_NAME} PROPERTIES
                         TIMEOUT 14400
                         PROCESSORS ${TEST_NP}
                         WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
                         LABELS "regression"
                         ATTACHED_FILES "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.log"
                         ATTACHED_FILES_ON_FAIL "${CURRENT_TEST_BINARY_DIR}/Backtrace.0")
endfunction(add_test_r)

# Regression tests excluded from CI
function(add_test_re TEST_NAME)
    add_test_r(${TEST_NAME})
    set_tests_properties(${TEST_NAME} PROPERTIES LABELS "regression;no_ci")
endfunction(add_test_re)

# Regression test and excluded from CI with dependency
function(add_test_red TEST_NAME TEST_DEPENDENCY)
    add_test_re(${TEST_NAME})
    set_tests_properties(${TEST_NAME} PROPERTIES FIXTURES_REQUIRED fixture_${TEST_DEPENDENCY})
    set_tests_properties(${TEST_DEPENDENCY} PROPERTIES FIXTURES_SETUP fixture_${TEST_DEPENDENCY})
endfunction(add_test_red)

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
      set(RUN_COMMAND_${GRID_SIZE} "${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${amr_wind_exe_name} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${GRID_SIZE}/${TEST_NAME}.inp")
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
endfunction(add_test_v)

# Standard unit test
function(add_test_u TEST_NAME)
    setup_test()
    set(TEST_NP 1)
    if(AMR_WIND_ENABLE_MPI)
      set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${TEST_NP} ${MPIEXEC_PREFLAGS}")
    else()
      unset(MPI_COMMANDS)
    endif()
    add_test(${TEST_NAME} bash -c "${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${amr_wind_unit_test_exe_name}")
    set_tests_properties(${TEST_NAME} PROPERTIES
                         TIMEOUT 500
                         PROCESSORS ${TEST_NP}
                         WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
                         LABELS "unit")
endfunction(add_test_u)

# Add test that also checks post-processing scripts: sampling
function(add_test_pps TEST_NAME FORMAT NINT NSTEPS)
    setup_test()
    set(ADDL_RUNTIME_OPTIONS "sampling.output_format=${FORMAT} sampling.output_interval=${NINT} time.plot_interval=-1 time.checkpoint_interval=-1 time.max_step=${NSTEPS}")
    add_test("${TEST_NAME}_sampling_${FORMAT}" bash -c "set -o pipefail && ${MPI_COMMANDS} ${CMAKE_BINARY_DIR}/${amr_wind_exe_name} ${MPIEXEC_POSTFLAGS} ${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}.inp ${RUNTIME_OPTIONS} ${ADDL_RUNTIME_OPTIONS} 2>&1 | tee ${TEST_NAME}_sampling_${FORMAT}.log && ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/sampling_${TEST_NAME}_${FORMAT}.py")
    # Set properties for test
    set_tests_properties("${TEST_NAME}_sampling_${FORMAT}" PROPERTIES
                          TIMEOUT 14400
                          PROCESSORS ${TEST_NP}
                          WORKING_DIRECTORY "${CURRENT_TEST_BINARY_DIR}/"
                          LABELS "post_processing;no_ci"
                          ATTACHED_FILES "${CURRENT_TEST_BINARY_DIR}/${TEST_NAME}_sampling_${FORMAT}.log;plot_sampling_${FORMAT}.pdf")
endfunction(add_test_pps)

#=============================================================================
# Unit tests
#=============================================================================
add_test_u(unit_tests)

#=============================================================================
# Regression tests
#=============================================================================
add_test_r(abl_godunov)
add_test_r(abl_mol)
add_test_r(boussinesq_bubble_godunov)
add_test_r(freestream_godunov)
add_test_r(tgv_godunov)

#=============================================================================
# Regression tests excluded from CI
#=============================================================================
add_test_re(abl_godunov_mpl)
add_test_re(abl_godunov_mpl_amr)
add_test_re(abl_godunov_cn)
add_test_re(abl_godunov_explicit)
add_test_re(abl_godunov_nolim)
add_test_re(abl_godunov_plm)
add_test_re(abl_godunov_noncons)
add_test_re(abl_bds)
add_test_re(abl_godunov_static_refinement)
add_test_re(abl_godunov_static_refinement_rr4)
add_test_re(abl_godunov_scalar_velocity_solve)
add_test_re(abl_godunov_segregated_velocity_solve)
add_test_re(abl_godunov_timetable)
add_test_re(abl_godunov_geostrophic_timetable)
add_test_re(abl_ksgsm84_godunov)
add_test_re(abl_mol_cn)
add_test_re(abl_mol_explicit)
add_test_re(abl_sampling)
add_test_re(abl_subvolume)
add_test_re(abl_stable)
add_test_re(abl_stable_different_roughness)
add_test_re(abl_unstable)
add_test_re(abl_unstable_constant_wall_model)
add_test_re(abl_unstable_local_wall_model)
add_test_re(abl_unstable_schumann_wall_model)
add_test_re(abl_surf_temp_timetable)
add_test_re(abl_anelastic)
add_test_re(abl_multiphase_laminar)
add_test_re(abl_multiphase_laminar_2planes)
add_test_re(abl_waves_terrain)
add_test_re(act_abl_joukowskydisk)
add_test_re(act_abl_uniformctdisk)
add_test_re(act_fixed_wing)
add_test_re(act_fixed_wing_fllc)
add_test_re(act_moving_wing)
add_test_re(act_pitching_wing_2D)
add_test_re(act_flat_plate)
add_test_re(boussinesq_bubble_mol)
add_test_re(uniform_ct_disk)
add_test_re(uniform_ct_disk_gaussian)
add_test_re(uniform_ct_disk_dynamic_adaptation)
add_test_re(joukowsky_disk)
add_test_re(joukowsky_disk_awc)
add_test_re(channel_kwsst)
add_test_re(channel_kwsst_sust)
add_test_re(channel_kwsstiddes)
add_test_re(channel_godunov_laminar)
add_test_re(channel_smagorinsky_analytical)
add_test_re(halfchannel_zerogradient)
add_test_re(halfchannel_symmetricwall)
add_test_re(ekman_spiral)
add_test_re(rayleigh_taylor_godunov)
add_test_re(rayleigh_taylor_mol)
add_test_re(tgv_godunov_plm)
add_test_re(tgv_mol)
add_test_re(vortex_patch_godunov)
add_test_re(zalesak_disk_godunov)
add_test_re(dam_break_godunov)
add_test_re(sloshing_tank)
add_test_re(abl_godunov_weno)
add_test_re(abl_godunov_wenoz_fixedpt)
add_test_re(abl_godunov_ppm)
add_test_re(abl_amd_wenoz)
add_test_re(ib_ctv_godunov_weno)
add_test_re(ib_cylinder_Re_300)
add_test_re(ib_sphere_Re_100)
add_test_re(vortex_ring_collision)
add_test_re(fat_cored_vortex_ring)
add_test_re(abl_bndry_output_native)
add_test_re(abl_bndry_output_amr_native)
add_test_re(vortex_patch_scalar_vel)
add_test_re(zalesak_disk_scalar_vel)
add_test_re(rain_drop)
add_test_re(inertial_drop)
add_test_re(ow_current_wing)
add_test_re(ow_current_bathymetry)
add_test_re(scalar_advection_uniform)
add_test_re(scalar_advection_refined)
add_test_re(freestream_bds)
add_test_re(inflow_godunov_amr)
add_test_re(inflow_bds_amr)
add_test_re(linear_godunov_amr)
add_test_re(linear_bds_amr)
add_test_re(hbl_godunov)
add_test_re(vortex_dipole_wall_collision)
add_test_re(vortex_dipole_inout_refine)
add_test_re(burggraf_flow)
add_test_re(abl_godunov_rayleigh_damping)
add_test_re(terrain_box)
add_test_re(terrain_box_amr)
add_test_re(forest_drag)
add_test_re(box_refinement)
add_test_re(cylinder_refinement)
add_test_re(freestream_godunov_inout)
add_test_re(freestream_godunov_inout_refine1)
add_test_re(freestream_godunov_inout_refine2)
add_test_re(ctv_godunov_plm)
add_test_re(ctv_bds)
add_test_re(ctv_mol_mesh_map)
add_test_re(ctv_mol_mesh_map_explicit)
add_test_re(abl_specifiedmol_neutral)
add_test_re(abl_specifiedmol_stable)
add_test_re(abl_specifiedmol_unstable)
add_test_re(udf_refinement)

if("${AMR_WIND_PRECISION}" STREQUAL "DOUBLE")
  add_test_re(rankine)
  add_test_re(rankine-sym)
  add_test_re(ow_linear)
  add_test_re(ow_linear_init_waves)
  add_test_re(ow_stokes)
endif()

if(AMR_WIND_ENABLE_NETCDF)
  add_test_re(abl_bndry_output)
  add_test_re(abl_bndry_output_amr_inflow)
  add_test_re(abl_bndry_output_amr_upper)
  add_test_re(abl_virtual_lidar)
  add_test_re(abl_meso_input_dpa)
  add_test_re(abl_meso_input_ipa)
  add_test_re(abl_meso_tendency)
  add_test_re(abl_kosovic_neutral)
  add_test_re(abl_sampling_netcdf)
  add_test_re(abl_kosovic_neutral_ib)
  add_test_re(nrel_precursor)
  add_test_re(abl_wallrans_neutral)
endif()

if(AMR_WIND_ENABLE_MASA)
  add_test_re(mms_godunov)
  add_test_re(mms_godunov_plm)
  add_test_re(mms_bds)
  add_test_re(mms_mol)
endif()

if(AMR_WIND_ENABLE_HYPRE)
  add_test_re(abl_godunov_hypre)
  add_test_re(channel_kwsst_hypre)
  add_test_re(channel_mol_mesh_map_x)
  add_test_re(channel_mol_mesh_map_y)
  add_test_re(channel_mol_mesh_map_z)
  add_test_re(channel_mol_mesh_map_x_seg_vel_solve)
endif()

if(AMR_WIND_ENABLE_HDF5)
  add_test_re(abl_stable_hdf5)
  if(AMR_WIND_ENABLE_HDF5_ZFP)
    add_test_re(abl_stable_zfp_hdf5)
  endif()
endif()

if(AMR_WIND_ENABLE_ASCENT)
  add_test_re(abl_godunov_ascent)
endif()

if(AMR_WIND_ENABLE_W2A)
  add_test_re(ow_w2a)
  add_test_re(ow_w2a_nwt_2d)
  add_test_re(ow_w2a_nwt_3d)
  add_test_re(ow_w2a_nonperiodic)
  add_test_re(abl_multiphase_w2a)
  add_test_red(abl_w2a_terrain abl_bndry_output_native)
endif()

if(AMR_WIND_ENABLE_OPENFAST AND AMR_WIND_OPENFAST_VERSION VERSION_LESS "4")
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
  if(AMR_WIND_ENABLE_NETCDF)
    add_test_red(act_uniform_alm_restart act_uniform_alm)
  endif()
endif()

if(AMR_WIND_ENABLE_KYNEMA)
  add_test_re(act_kynema_alm)
endif()

if(AMR_WIND_ENABLE_FFT)
  add_test_re(abl_godunov_nofft)
endif()
#=============================================================================
# Regression tests excluded from CI with a test dependency
#=============================================================================
if(AMR_WIND_TEST_WITH_FCOMPARE)
  add_test_red(abl_bndry_input_native abl_bndry_output_native)
  add_test_red(abl_bndry_input_native_inout abl_bndry_output_native)
  add_test_red(abl_godunov_restart abl_godunov)
  add_test_red(abl_bndry_input_amr_native abl_bndry_output_native)
  add_test_red(abl_bndry_input_amr_native_xhi abl_bndry_output_native)
  add_test_red(abl_bndry_input_amr_native_mlbc abl_bndry_output_amr_native)
  add_test_red(abl_godunov_forcetimetable abl_godunov_timetable)
  add_test_red(abl_multiphase_laminar_input abl_multiphase_laminar)
  add_test_red(abl_multiphase_laminar_inout abl_multiphase_laminar_2planes)
endif()

if(AMR_WIND_ENABLE_NETCDF)
  add_test_red(abl_bndry_input abl_bndry_output)
  add_test_red(abl_bndry_input_init abl_bndry_output)
  add_test_red(abl_bndry_input_amr abl_bndry_output)
  add_test_red(abl_bndry_input_amr_inflow abl_bndry_output_amr_inflow)
  add_test_red(abl_bndry_input_amr_upper abl_bndry_output_amr_upper)
  add_test_red(nrel_terrain nrel_precursor)
  add_test_red(nrel_terrain_amr nrel_precursor)
endif()

#=============================================================================
# Verification tests
#=============================================================================
if(AMR_WIND_ENABLE_MASA)
  set(LIST_OF_GRID_SIZES 8 16 32 64)
  add_test_v(mms "${LIST_OF_GRID_SIZES}")
endif()

#=============================================================================
# Postprocessing tests
#=============================================================================
if (AMR_WIND_TEST_WITH_PYTHON)
  add_test_pps(dam_break_godunov native 30 150)
  add_test_pps(dam_break_godunov ascii 30 150)
  if(AMR_WIND_ENABLE_NETCDF)
    add_test_pps(dam_break_godunov netcdf 30 150)
  endif()
endif()
