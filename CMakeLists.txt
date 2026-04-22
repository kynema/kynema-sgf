############################ BASE ######################################

cmake_minimum_required (VERSION 3.20 FATAL_ERROR)
project(AMR-Wind CXX C)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(amr-wind-utils)

########################## OPTIONS #####################################

#General options for the project
option(AMR_WIND_ENABLE_CLANG_TIDY "Compile with clang-tidy static analysis" OFF)
option(AMR_WIND_ENABLE_CPPCHECK "Enable cppcheck static analysis target" OFF)
option(AMR_WIND_ENABLE_FCOMPARE "Enable building fcompare when not testing" OFF)
option(AMR_WIND_ENABLE_DOCUMENTATION "Enable documentation target" OFF)

#Enabling tests overrides the executable options
option(AMR_WIND_ENABLE_UNIT_TESTS "Enable unit testing" ON)
option(AMR_WIND_ENABLE_TESTS "Enable testing suite" OFF)
option(AMR_WIND_TEST_WITH_FCOMPARE "Check test plots against gold files" OFF)
option(AMR_WIND_SAVE_GOLDS "Provide a directory in which to save golds during testing" OFF)
option(AMR_WIND_ENABLE_FPE_TRAP_FOR_TESTS "Enable FPE trapping in tests" ON)
option(AMR_WIND_TEST_WITH_PYTHON "Enable testing relying on python, like postprocessing" OFF)

#Use internal AMReX and AMReX-Hydro by default
option(AMR_WIND_USE_INTERNAL_AMREX "Use the AMReX submodule" ON)
option(AMR_WIND_USE_INTERNAL_AMREX_HYDRO "Use the AMReX-Hydro submodule" ON)

#Options for the executable
option(AMR_WIND_ENABLE_MPI "Enable MPI" OFF)
option(AMR_WIND_ENABLE_OPENMP "Enable OpenMP" OFF)
option(AMR_WIND_ENABLE_CUDA "Enable CUDA" OFF)
option(AMR_WIND_ENABLE_CUDA_RDC "Enable CUDA relocatable device code" OFF)
option(AMR_WIND_ENABLE_ROCM "Enable ROCm/HIP" OFF)
option(AMR_WIND_ENABLE_SYCL "Enable Intel OneAPI SyCL" OFF)
option(AMR_WIND_ENABLE_TINY_PROFILE "Enable AMReX TinyProfile support" OFF)
set(AMR_WIND_PRECISION "DOUBLE" CACHE STRING "Floating point precision SINGLE or DOUBLE")

# Third party libraries
option(AMR_WIND_ENABLE_HDF5 "Enable HDF5 library" OFF)
option(AMR_WIND_ENABLE_HDF5_ZFP "Enable ZFP compression in HDF5 library" OFF)
option(AMR_WIND_ENABLE_NETCDF "Enable NetCDF library" OFF)
option(AMR_WIND_ENABLE_MASA "Enable MASA library" OFF)
option(AMR_WIND_ENABLE_HELICS "Enable HELICS library" OFF)
option(AMR_WIND_ENABLE_HYPRE "Enable HYPRE integration" OFF)
option(AMR_WIND_ENABLE_OPENFAST "Enable OpenFAST integration" OFF)
option(AMR_WIND_ENABLE_KYNEMA "Enable Kynema integration" OFF)
option(AMR_WIND_ENABLE_ASCENT "Enable Ascent visualization library" OFF)
option(AMR_WIND_ENABLE_UMPIRE "Enable Umpire GPU memory pools" OFF)
option(AMR_WIND_ENABLE_W2A "Enable Waves2AMR library" OFF)
option(AMR_WIND_ENABLE_FFT "Enable FFT solver for MAC projection" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if((AMR_WIND_ENABLE_SYCL OR AMR_WIND_ENABLE_ROCM OR AMR_WIND_ENABLE_CUDA) AND AMR_WIND_ENABLE_OPENMP)
  message(FATAL_ERROR "Using OpenMP with a GPU backend is unsupported.")
endif()

if(AMR_WIND_ENABLE_CUDA)
  enable_language(CUDA)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "12.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION} which is unsupported."
      "Please use CUDA toolkit version 12.0 or newer.")
  endif()
endif()

if(AMR_WIND_ENABLE_ROCM)
  find_package(HIP REQUIRED)
  if(NOT DEFINED AMD_ARCH)
    # Set default AMD architectures (based on Frontier)
    set(AMD_ARCH "gfx90a")
  endif()
endif()

if(AMR_WIND_TEST_WITH_FCOMPARE)
  set(AMR_WIND_ENABLE_FCOMPARE ON)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  set(AMR_WIND_ENABLE_FPE_TRAP_FOR_TESTS OFF)
  message(WARNING "Disabling FPE trapping for tests when using AppleClang.")
endif()

if(AMR_WIND_ENABLE_KYNEMA AND (NOT AMR_WIND_ENABLE_MPI))
  message(FATAL_ERROR "Enabling Kyenma currently requires MPI to be enabled.")
endif()

if(AMR_WIND_ENABLE_OPENFAST)
  if(NOT AMR_WIND_ENABLE_MPI)
    message(FATAL_ERROR "Enabling OpenFAST currently requires MPI to be enabled.")
  endif()
  enable_language(Fortran)
endif()

########################### AMReX #####################################
init_amrex()
init_amrex_hydro()

########################## Waves2AMR ##################################
if(AMR_WIND_ENABLE_W2A)
  init_waves2amr()
endif()

########################### AMR-Wind #####################################
message(STATUS "AMR-Wind Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

include(set_rpath)

#Create target names
set(amr_wind_lib_name "amrwind_obj")
set(amr_wind_exe_name "amr_wind")
set(amr_wind_unit_test_exe_name "${amr_wind_exe_name}_unit_tests")
set(aw_api_lib "amrwind_api")

add_library(${amr_wind_lib_name} OBJECT)
add_library(${aw_api_lib})
add_executable(${amr_wind_exe_name})

init_code_checks()
if(CLANG_TIDY_EXE)
  set_target_properties(${amr_wind_lib_name} ${aw_api_lib} ${amr_wind_exe_name}
                        PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
endif()

include(set_compile_flags)

if(AMR_WIND_ENABLE_W2A)
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_W2A)
  target_link_libraries(${amr_wind_lib_name} PUBLIC Waves2AMR::Waves2AMR)
endif()

if(AMR_WIND_ENABLE_NETCDF)
  set(CMAKE_PREFIX_PATH ${NETCDF_DIR} ${CMAKE_PREFIX_PATH})
  find_package(NetCDF REQUIRED)
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_NETCDF)
  target_link_libraries_system(${amr_wind_lib_name} PUBLIC NetCDF::NetCDF)
endif()

if(AMR_WIND_ENABLE_UMPIRE)
  set(CMAKE_PREFIX_PATH ${UMPIRE_DIR} ${CMAKE_PREFIX_PATH})
  find_package(Umpire REQUIRED)
  target_include_directories(${amr_wind_lib_name} SYSTEM ${UMPIRE_INCLUDE_DIRS})
endif()

if(AMR_WIND_ENABLE_HDF5)
  set(CMAKE_PREFIX_PATH ${HDF5_DIR} ${CMAKE_PREFIX_PATH})
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_HDF5)
endif()

if(AMR_WIND_ENABLE_HDF5_ZFP)
  if(NOT AMR_WIND_ENABLE_HDF5)
    message(FATAL_ERROR "HDF5 must be enabled to use ZFP compression in HDF5 plot files")
  endif()
  set(CMAKE_PREFIX_PATH ${H5Z_ZFP_DIR} ${CMAKE_PREFIX_PATH})
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_HDF5_ZFP)
endif()

if(AMR_WIND_ENABLE_MASA)
  set(CMAKE_PREFIX_PATH ${MASA_DIR} ${CMAKE_PREFIX_PATH})
  find_package(MASA REQUIRED)
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_MASA)
  target_link_libraries(${amr_wind_lib_name} PUBLIC MASA::MASA)
endif()

if(AMR_WIND_ENABLE_OPENFAST)
  set(CMAKE_PREFIX_PATH ${OPENFAST_DIR} ${CMAKE_PREFIX_PATH})
  find_package(OpenFAST REQUIRED)
  include(CheckIncludeFileCXX)
  set(CMAKE_REQUIRED_INCLUDES "${OpenFAST_INCLUDE_DIRS}")
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_OPENFAST)
  target_include_directories(${amr_wind_lib_name} PUBLIC ${OpenFAST_INCLUDE_DIRS})
  target_link_libraries(${amr_wind_lib_name} PUBLIC ${OpenFAST_LIBRARIES})

  set(AMR_WIND_OPENFAST_VERSION "3" CACHE STRING "OpenFAST major version")
  if((AMR_WIND_OPENFAST_VERSION STREQUAL "master") OR (AMR_WIND_OPENFAST_VERSION STREQUAL "main") OR (AMR_WIND_OPENFAST_VERSION STREQUAL "develop"))
    set(AMR_WIND_OPENFAST_VERSION "4.1")
  endif()
  if(NOT AMR_WIND_OPENFAST_VERSION VERSION_LESS "5")
    message(FATAL_ERROR "AMR_WIND_OPENFAST_VERSION must be less than 5.")
  endif()
  if(AMR_WIND_OPENFAST_VERSION VERSION_GREATER_EQUAL "4" AND AMR_WIND_OPENFAST_VERSION VERSION_LESS "4.1")
    message(FATAL_ERROR "AMR-Wind does not support OpenFAST version 4.0. Either update to 4.1 or revert to 3.")
  endif()
  message(STATUS "AMR_WIND_OPENFAST_VERSION = ${AMR_WIND_OPENFAST_VERSION}")
  string(REPLACE "." ";" OPENFAST_VERSION_LIST ${AMR_WIND_OPENFAST_VERSION})
  list(LENGTH OPENFAST_VERSION_LIST OPENFAST_VERSION_LIST_LENGTH)
  if(OPENFAST_VERSION_LIST_LENGTH GREATER_EQUAL 1)
    list(GET OPENFAST_VERSION_LIST 0 OPENFAST_VERSION_MAJOR)
  endif()
  if(OPENFAST_VERSION_LIST_LENGTH GREATER_EQUAL 2)
    list(GET OPENFAST_VERSION_LIST 1 OPENFAST_VERSION_MINOR)
  endif()
  target_compile_definitions(${amr_wind_lib_name} PUBLIC OPENFAST_VERSION_MAJOR=${OPENFAST_VERSION_MAJOR})
  target_compile_definitions(${amr_wind_lib_name} PUBLIC OPENFAST_VERSION_MINOR=${OPENFAST_VERSION_MINOR})
endif()

if(AMR_WIND_ENABLE_KYNEMA)
  set(CMAKE_PREFIX_PATH ${KYNEMA_DIR} ${CMAKE_PREFIX_PATH})
  find_package(Kynema REQUIRED)
  target_link_libraries(${amr_wind_lib_name} PUBLIC Kynema::kynema_library)
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_KYNEMA)
endif()

if(AMR_WIND_ENABLE_ASCENT)
  set(CMAKE_PREFIX_PATH ${ASCENT_DIR} ${CMAKE_PREFIX_PATH})
  find_package(Ascent REQUIRED)
  if(AMR_WIND_ENABLE_MPI)
    target_link_libraries(${amr_wind_lib_name} PUBLIC ascent::ascent_mpi)
  else()
    target_link_libraries(${amr_wind_lib_name} PUBLIC ascent::ascent)
  endif()
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_ASCENT)
endif()

if(AMR_WIND_ENABLE_HELICS)
  set(CMAKE_PREFIX_PATH ${HELICS_DIR} ${CMAKE_PREFIX_PATH})
  find_package(HELICS 3 REQUIRED)
  target_include_directories(${amr_wind_lib_name} PUBLIC ${HELICS_INSTALL_DIR}/include)
  target_link_libraries(${amr_wind_lib_name} PUBLIC HELICS::helicsCpp98)
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_HELICS)
endif()

if(AMR_WIND_ENABLE_FFT)
  target_compile_definitions(${amr_wind_lib_name} PUBLIC AMR_WIND_USE_FFT)
endif()

if(AMR_WIND_TEST_WITH_PYTHON)
  set(CMAKE_PREFIX_PATH ${PYTHON_DIR} ${CMAKE_PREFIX_PATH})
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
endif()

generate_version_info()
add_subdirectory(src)

target_link_libraries_system(${amr_wind_lib_name} PUBLIC AMReX::amrex AMReX-Hydro::amrex_hydro_api)
target_link_libraries(${amr_wind_exe_name} PRIVATE ${amr_wind_lib_name} AMReX-Hydro::amrex_hydro_api)
target_link_libraries(${aw_api_lib} PUBLIC ${amr_wind_lib_name} AMReX-Hydro::amrex_hydro_api)

set_target_properties(${amr_wind_lib_name} buildInfo${amr_wind_lib_name}
	  PROPERTIES POSITION_INDEPENDENT_CODE ON)

if(AMR_WIND_ENABLE_CUDA)
  set(awtargets "${amr_wind_lib_name};${amr_wind_exe_name}")
  foreach(tgt IN LISTS awtargets)
    set_cuda_build_properties(${tgt})
  endforeach()
endif()

if(AMR_WIND_ENABLE_UNIT_TESTS OR AMR_WIND_ENABLE_TESTS)
  add_executable(${amr_wind_unit_test_exe_name})
  if(CLANG_TIDY_EXE)
    set_target_properties(${amr_wind_unit_test_exe_name}
                          PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
  endif()
  add_subdirectory("submods/googletest")
  add_subdirectory("unit_tests")
  set_cuda_build_properties(${amr_wind_unit_test_exe_name})
endif()

add_subdirectory(tools)

if(AMR_WIND_ENABLE_TESTS)
  enable_testing()
  include(CTest)
  add_subdirectory(test)
endif()

if(AMR_WIND_ENABLE_DOCUMENTATION)
  add_subdirectory(docs)
endif()

#Define what we want to be installed during a make install
install(
  TARGETS ${amr_wind_exe_name}
  ${aw_api_lib} ${amr_wind_lib_name} buildInfo${amr_wind_lib_name}
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

if(AMR_WIND_ENABLE_UNIT_TESTS OR AMR_WIND_ENABLE_TESTS)
  install(
    TARGETS ${amr_wind_unit_test_exe_name}
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
