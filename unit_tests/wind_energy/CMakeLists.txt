target_sources(${kynema_sgf_unit_test_exe_name} PRIVATE
  # testing utilities
  abl_test_utils.cpp

  # test cases
  test_abl_init.cpp
  test_abl_src.cpp
  test_abl_stats.cpp
  test_abl_bc.cpp
  test_abl_src_timetable.cpp
  test_abl_terrain.cpp
  test_abl_forest.cpp
  )

if (AMR_WIND_ENABLE_NETCDF)
  target_sources(${kynema_sgf_unit_test_exe_name} PRIVATE
    test_abl_init_ncf.cpp
    )
endif()

add_subdirectory(actuator)
