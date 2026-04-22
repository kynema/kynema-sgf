target_sources(${amr_wind_lib_name} PRIVATE
  turbine_utils.cpp
  )

add_subdirectory(fast)
add_subdirectory(kynema)