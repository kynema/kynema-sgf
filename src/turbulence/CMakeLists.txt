target_sources(${kynema_sgf_lib_name} PRIVATE
  LaminarModel.cpp
  turb_utils.cpp
  )

add_subdirectory(LES)
add_subdirectory(RANS)
