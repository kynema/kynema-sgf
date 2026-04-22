target_sources(${kynema_sgf_unit_test_exe_name} PRIVATE
  test_refinement.cpp
  test_plane_averaging.cpp
  test_field_plane_averaging.cpp
  test_field_plane_averaging_fine.cpp
  test_second_moment.cpp
  test_sampling.cpp
  test_linear_interpolation.cpp
  test_integrals.cpp
  test_free_surface.cpp
  test_wave_energy.cpp
  test_diagnostics.cpp
  test_multilevelvector.cpp
  test_field_norms.cpp
  test_tensor_ops.cpp
  test_post_processing_time.cpp
  test_time_averaging.cpp
  )

if (AMR_WIND_ENABLE_NETCDF)
  target_sources(${kynema_sgf_unit_test_exe_name} PRIVATE
    test_ncutils.cpp
    )
endif()
