target_sources(${kynema_sgf_lib_name} PRIVATE
  PDEBase.cpp
  DiffusionOps.cpp
  )

add_subdirectory(icns)
add_subdirectory(temperature)
add_subdirectory(density)
add_subdirectory(tke)
add_subdirectory(sdr)
add_subdirectory(levelset)
add_subdirectory(vof)
add_subdirectory(passive_scalar)
