target_sources(${amr_wind_lib_name} PRIVATE
  LaminarModel.cpp
  turb_utils.cpp
  )

add_subdirectory(LES)
add_subdirectory(RANS)
