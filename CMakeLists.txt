############################ BASE ######################################

cmake_minimum_required (VERSION 3.20 FATAL_ERROR)
project(Kynema-SGF CXX C)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(kynema-sgf-utils)

########################## OPTIONS #####################################

#General options for the project
option(KYNEMA_SGF_ENABLE_CLANG_TIDY "Compile with clang-tidy static analysis" OFF)
option(KYNEMA_SGF_ENABLE_CPPCHECK "Enable cppcheck static analysis target" OFF)
option(KYNEMA_SGF_ENABLE_FCOMPARE "Enable building fcompare when not testing" OFF)
option(KYNEMA_SGF_ENABLE_DOCUMENTATION "Enable documentation target" OFF)

#Enabling tests overrides the executable options
option(KYNEMA_SGF_ENABLE_UNIT_TESTS "Enable unit testing" ON)
option(KYNEMA_SGF_ENABLE_TESTS "Enable testing suite" OFF)
option(KYNEMA_SGF_TEST_WITH_FCOMPARE "Check test plots against gold files" OFF)
option(KYNEMA_SGF_SAVE_GOLDS "Provide a directory in which to save golds during testing" OFF)
option(KYNEMA_SGF_ENABLE_FPE_TRAP_FOR_TESTS "Enable FPE trapping in tests" ON)
option(KYNEMA_SGF_TEST_WITH_PYTHON "Enable testing relying on python, like postprocessing" OFF)

#Use internal AMReX and AMReX-Hydro by default
option(KYNEMA_SGF_USE_INTERNAL_AMREX "Use the AMReX submodule" ON)
option(KYNEMA_SGF_USE_INTERNAL_AMREX_HYDRO "Use the AMReX-Hydro submodule" ON)

#Options for the executable
option(KYNEMA_SGF_ENABLE_MPI "Enable MPI" OFF)
option(KYNEMA_SGF_ENABLE_OPENMP "Enable OpenMP" OFF)
option(KYNEMA_SGF_ENABLE_CUDA "Enable CUDA" OFF)
option(KYNEMA_SGF_ENABLE_CUDA_RDC "Enable CUDA relocatable device code" OFF)
option(KYNEMA_SGF_ENABLE_ROCM "Enable ROCm/HIP" OFF)
option(KYNEMA_SGF_ENABLE_SYCL "Enable Intel OneAPI SyCL" OFF)
option(KYNEMA_SGF_ENABLE_TINY_PROFILE "Enable AMReX TinyProfile support" OFF)
set(KYNEMA_SGF_PRECISION "DOUBLE" CACHE STRING "Floating point precision SINGLE or DOUBLE")

# Third party libraries
option(KYNEMA_SGF_ENABLE_HDF5 "Enable HDF5 library" OFF)
option(KYNEMA_SGF_ENABLE_HDF5_ZFP "Enable ZFP compression in HDF5 library" OFF)
option(KYNEMA_SGF_ENABLE_NETCDF "Enable NetCDF library" OFF)
option(KYNEMA_SGF_ENABLE_MASA "Enable MASA library" OFF)
option(KYNEMA_SGF_ENABLE_HELICS "Enable HELICS library" OFF)
option(KYNEMA_SGF_ENABLE_HYPRE "Enable HYPRE integration" OFF)
option(KYNEMA_SGF_ENABLE_OPENFAST "Enable OpenFAST integration" OFF)
option(KYNEMA_SGF_ENABLE_KYNEMA "Enable Kynema integration" OFF)
option(KYNEMA_SGF_ENABLE_ASCENT "Enable Ascent visualization library" OFF)
option(KYNEMA_SGF_ENABLE_UMPIRE "Enable Umpire GPU memory pools" OFF)
option(KYNEMA_SGF_ENABLE_W2A "Enable Waves2AMR library" OFF)
option(KYNEMA_SGF_ENABLE_FFT "Enable FFT solver for MAC projection" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if((KYNEMA_SGF_ENABLE_SYCL OR KYNEMA_SGF_ENABLE_ROCM OR KYNEMA_SGF_ENABLE_CUDA) AND KYNEMA_SGF_ENABLE_OPENMP)
  message(FATAL_ERROR "Using OpenMP with a GPU backend is unsupported.")
endif()

if(KYNEMA_SGF_ENABLE_CUDA)
  enable_language(CUDA)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "12.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION} which is unsupported."
      "Please use CUDA toolkit version 12.0 or newer.")
  endif()
endif()

if(KYNEMA_SGF_ENABLE_ROCM)
  find_package(HIP REQUIRED)
  if(NOT DEFINED AMD_ARCH)
    # Set default AMD architectures (based on Frontier)
    set(AMD_ARCH "gfx90a")
  endif()
endif()

if(KYNEMA_SGF_TEST_WITH_FCOMPARE)
  set(KYNEMA_SGF_ENABLE_FCOMPARE ON)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  set(KYNEMA_SGF_ENABLE_FPE_TRAP_FOR_TESTS OFF)
  message(WARNING "Disabling FPE trapping for tests when using AppleClang.")
endif()

if(KYNEMA_SGF_ENABLE_KYNEMA AND (NOT KYNEMA_SGF_ENABLE_MPI))
  message(FATAL_ERROR "Enabling Kyenma currently requires MPI to be enabled.")
endif()

if(KYNEMA_SGF_ENABLE_OPENFAST)
  if(NOT KYNEMA_SGF_ENABLE_MPI)
    message(FATAL_ERROR "Enabling OpenFAST currently requires MPI to be enabled.")
  endif()
  enable_language(Fortran)
endif()

########################### AMReX #####################################
init_amrex()
init_amrex_hydro()

########################## Waves2AMR ##################################
if(KYNEMA_SGF_ENABLE_W2A)
  init_waves2amr()
endif()

########################### Kynema-SGF #####################################
message(STATUS "Kynema-SGF Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

include(set_rpath)

#Create target names
set(kynema_sgf_lib_name "kynema_sgf_obj")
set(kynema_sgf_exe_name "kynema_sgf")
set(kynema_sgf_unit_test_exe_name "${kynema_sgf_exe_name}_unit_tests")
set(ks_api_lib "kynema_sgf_api")

add_library(${kynema_sgf_lib_name} OBJECT)
add_library(${ks_api_lib})
add_executable(${kynema_sgf_exe_name})

init_code_checks()
if(CLANG_TIDY_EXE)
  set_target_properties(${kynema_sgf_lib_name} ${ks_api_lib} ${kynema_sgf_exe_name}
                        PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
endif()

include(set_compile_flags)

if(KYNEMA_SGF_ENABLE_W2A)
  target_compile_definitions(${kynema_sgf_lib_name} PUBLIC KYNEMA_SGF_USE_W2A)
  target_link_libraries(${kynema_sgf_lib_name} PUBLIC Waves2AMR::Waves2AMR)
endif()

if(KYNEMA_SGF_ENABLE_NETCDF)
  set(CMAKE_PREFIX_PATH ${NETCDF_DIR} ${CMAKE_PREFIX_PATH})
  find_package(NetCDF REQUIRED)
  target_compile_definitions(${kynema_sgf_lib_name} PUBLIC KYNEMA_SGF_USE_NETCDF)
  target_link_libraries_system(${kynema_sgf_lib_name} PUBLIC NetCDF::NetCDF)
endif()

if(KYNEMA_SGF_ENABLE_UMPIRE)
  set(CMAKE_PREFIX_PATH ${UMPIRE_DIR} ${CMAKE_PREFIX_PATH})
  find_package(Umpire REQUIRED)
  target_include_directories(${kynema_sgf_lib_name} SYSTEM ${UMPIRE_INCLUDE_DIRS})
endif()

if(KYNEMA_SGF_ENABLE_HDF5)
  set(CMAKE_PREFIX_PATH ${HDF5_DIR} ${CMAKE_PREFIX_PATH})
  target_compile_definitions(${kynema_sgf_lib_name} PUBLIC KYNEMA_SGF_USE_HDF5)
endif()

if(KYNEMA_SGF_ENABLE_HDF5_ZFP)
  if(NOT KYNEMA_SGF_ENABLE_HDF5)
    message(FATAL_ERROR "HDF5 must be enabled to use ZFP compression in HDF5 plot files")
  endif()
  set(CMAKE_PREFIX_PATH ${H5Z_ZFP_DIR} ${CMAKE_PREFIX_PATH})
  target_compile_definitions(${kynema_sgf_lib_name} PUBLIC KYNEMA_SGF_USE_HDF5_ZFP)
endif()

if(KYNEMA_SGF_ENABLE_MASA)
  set(CMAKE_PREFIX_PATH ${MASA_DIR} ${CMAKE_PREFIX_PATH})
  find_package(MASA REQUIRED)
  target_compile_definitions(${kynema_sgf_lib_name} PUBLIC KYNEMA_SGF_USE_MASA)
  target_link_libraries(${kynema_sgf_lib_name} PUBLIC MASA::MASA)
endif()

if(KYNEMA_SGF_ENABLE_OPENFAST)
  set(CMAKE_PREFIX_PATH ${OPENFAST_DIR} ${CMAKE_PREFIX_PATH})
  find_package(OpenFAST REQUIRED)
  include(CheckIncludeFileCXX)
  set(CMAKE_REQUIRED_INCLUDES "${OpenFAST_INCLUDE_DIRS}")
  target_compile_definitions(${kynema_sgf_lib_name} PUBLIC KYNEMA_SGF_USE_OPENFAST)
  target_include_directories(${kynema_sgf_lib_name} PUBLIC ${OpenFAST_INCLUDE_DIRS})
  target_link_libraries(${kynema_sgf_lib_name} PUBLIC ${OpenFAST_LIBRARIES})

  set(KYNEMA_SGF_OPENFAST_VERSION "3" CACHE STRING "OpenFAST major version")
  if((KYNEMA_SGF_OPENFAST_VERSION STREQUAL "master") OR (KYNEMA_SGF_OPENFAST_VERSION STREQUAL "main") OR (KYNEMA_SGF_OPENFAST_VERSION STREQUAL "develop"))
    set(KYNEMA_SGF_OPENFAST_VERSION "4.1")
  endif()
  if(NOT KYNEMA_SGF_OPENFAST_VERSION VERSION_LESS "5")
    message(FATAL_ERROR "KYNEMA_SGF_OPENFAST_VERSION must be less than 5.")
  endif()
  if(KYNEMA_SGF_OPENFAST_VERSION VERSION_GREATER_EQUAL "4" AND KYNEMA_SGF_OPENFAST_VERSION VERSION_LESS "4.1")
    message(FATAL_ERROR "Kynema-SGF does not support OpenFAST version 4.0. Either update to 4.1 or revert to 3.")
  endif()
  message(STATUS "KYNEMA_SGF_OPENFAST_VERSION = ${KYNEMA_SGF_OPENFAST_VERSION}")
  string(REPLACE "." ";" OPENFAST_VERSION_LIST ${KYNEMA_SGF_OPENFAST_VERSION})
  list(LENGTH OPENFAST_VERSION_LIST OPENFAST_VERSION_LIST_LENGTH)
  if(OPENFAST_VERSION_LIST_LENGTH GREATER_EQUAL 1)
    list(GET OPENFAST_VERSION_LIST 0 OPENFAST_VERSION_MAJOR)
  endif()
  if(OPENFAST_VERSION_LIST_LENGTH GREATER_EQUAL 2)
    list(GET OPENFAST_VERSION_LIST 1 OPENFAST_VERSION_MINOR)
  endif()
  target_compile_definitions(${kynema_sgf_lib_name} PUBLIC OPENFAST_VERSION_MAJOR=${OPENFAST_VERSION_MAJOR})
  target_compile_definitions(${kynema_sgf_lib_name} PUBLIC OPENFAST_VERSION_MINOR=${OPENFAST_VERSION_MINOR})
endif()

if(KYNEMA_SGF_ENABLE_KYNEMA)
  set(CMAKE_PREFIX_PATH ${KYNEMA_DIR} ${CMAKE_PREFIX_PATH})
  find_package(Kynema REQUIRED)
  target_link_libraries(${kynema_sgf_lib_name} PUBLIC Kynema::kynema_library)
  target_compile_definitions(${kynema_sgf_lib_name} PUBLIC KYNEMA_SGF_USE_KYNEMA)
endif()

if(KYNEMA_SGF_ENABLE_ASCENT)
  set(CMAKE_PREFIX_PATH ${ASCENT_DIR} ${CMAKE_PREFIX_PATH})
  find_package(Ascent REQUIRED)
  if(KYNEMA_SGF_ENABLE_MPI)
    target_link_libraries(${kynema_sgf_lib_name} PUBLIC ascent::ascent_mpi)
  else()
    target_link_libraries(${kynema_sgf_lib_name} PUBLIC ascent::ascent)
  endif()
  target_compile_definitions(${kynema_sgf_lib_name} PUBLIC KYNEMA_SGF_USE_ASCENT)
endif()

if(KYNEMA_SGF_ENABLE_HELICS)
  set(CMAKE_PREFIX_PATH ${HELICS_DIR} ${CMAKE_PREFIX_PATH})
  find_package(HELICS 3 REQUIRED)
  target_include_directories(${kynema_sgf_lib_name} PUBLIC ${HELICS_INSTALL_DIR}/include)
  target_link_libraries(${kynema_sgf_lib_name} PUBLIC HELICS::helicsCpp98)
  target_compile_definitions(${kynema_sgf_lib_name} PUBLIC KYNEMA_SGF_USE_HELICS)
endif()

if(KYNEMA_SGF_ENABLE_FFT)
  target_compile_definitions(${kynema_sgf_lib_name} PUBLIC KYNEMA_SGF_USE_FFT)
endif()

if(KYNEMA_SGF_TEST_WITH_PYTHON)
  set(CMAKE_PREFIX_PATH ${PYTHON_DIR} ${CMAKE_PREFIX_PATH})
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
endif()

generate_version_info()
add_subdirectory(src)

target_link_libraries_system(${kynema_sgf_lib_name} PUBLIC AMReX::amrex AMReX-Hydro::amrex_hydro_api)
target_link_libraries(${kynema_sgf_exe_name} PRIVATE ${kynema_sgf_lib_name} AMReX-Hydro::amrex_hydro_api)
target_link_libraries(${ks_api_lib} PUBLIC ${kynema_sgf_lib_name} AMReX-Hydro::amrex_hydro_api)

set_target_properties(${kynema_sgf_lib_name} buildInfo${kynema_sgf_lib_name}
	  PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(KYNEMA_SGF_ENABLE_CUDA)
  set(kstargets "${kynema_sgf_lib_name};${kynema_sgf_exe_name}")
  foreach(tgt IN LISTS kstargets)
    set_cuda_build_properties(${tgt})
  endforeach()
endif()

if(KYNEMA_SGF_ENABLE_UNIT_TESTS OR KYNEMA_SGF_ENABLE_TESTS)
  add_executable(${kynema_sgf_unit_test_exe_name})
  if(CLANG_TIDY_EXE)
    set_target_properties(${kynema_sgf_unit_test_exe_name}
                          PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
  endif()
  add_subdirectory("submods/googletest")
  add_subdirectory("unit_tests")
  set_cuda_build_properties(${kynema_sgf_unit_test_exe_name})
endif()

add_subdirectory(tools)

if(KYNEMA_SGF_ENABLE_TESTS)
  enable_testing()
  include(CTest)
  add_subdirectory(test)
endif()

if(KYNEMA_SGF_ENABLE_DOCUMENTATION)
  add_subdirectory(docs)
endif()

#Define what we want to be installed during a make install
install(
  TARGETS ${kynema_sgf_exe_name}
  ${ks_api_lib} ${kynema_sgf_lib_name} buildInfo${kynema_sgf_lib_name}
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(KYNEMA_SGF_ENABLE_UNIT_TESTS OR KYNEMA_SGF_ENABLE_TESTS)
  install(
    TARGETS ${kynema_sgf_unit_test_exe_name}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

install(
  DIRECTORY ${PROJECT_SOURCE_DIR}/src
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.H")

install(
  EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

configure_package_config_file(
  cmake/${PROJECT_NAME}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})
install(FILES
  ${PROJECT_SOURCE_DIR}/cmake/FindMASA.cmake
  ${PROJECT_SOURCE_DIR}/cmake/FindNetCDF.cmake
  DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/Modules)
