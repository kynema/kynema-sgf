target_sources(${amr_wind_lib_name}
  PRIVATE
  incflo_advance.cpp
  incflo.cpp
  incflo_compute_dt.cpp
  incflo_regrid.cpp
  helics.cpp
  CFDSim.cpp
  )

target_sources(${amr_wind_exe_name} PRIVATE main.cpp)

target_include_directories(${amr_wind_lib_name} PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

add_subdirectory(core)
add_subdirectory(boundary_conditions)
add_subdirectory(convection)
add_subdirectory(diffusion)
add_subdirectory(projection)
add_subdirectory(setup)
add_subdirectory(utilities)
add_subdirectory(wind_energy)
add_subdirectory(equation_systems)
add_subdirectory(transport_models)
add_subdirectory(turbulence)
add_subdirectory(physics)
add_subdirectory(overset)
add_subdirectory(immersed_boundary)
add_subdirectory(mesh_mapping_models)
add_subdirectory(ocean_waves)
add_subdirectory(eos_models)

include(AMReXBuildInfo)
generate_buildinfo(${amr_wind_lib_name} ${PROJECT_SOURCE_DIR})

configure_file("${PROJECT_SOURCE_DIR}/cmake/AMRWindVersion.H.in"
  "${CMAKE_CURRENT_BINARY_DIR}/AMRWindVersion.H" @ONLY)

file(GLOB_RECURSE amr_wind_header_list ${CMAKE_CURRENT_SOURCE_DIR} "*.H")
add_custom_target(amr_wind_headers SOURCES ${amr_wind_header_list})

install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/AMRWindVersion.H
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
