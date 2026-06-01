target_sources(${kynema_sgf_lib_name}
   PRIVATE
      #C++
      BCInterface.cpp
      FixedGradientBC.cpp
      scalar_bcs.cpp
      MassInflowOutflowBC.cpp
   )

add_subdirectory(wall_models)
add_subdirectory(field_boundary_fill)