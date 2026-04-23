target_sources(
  ${kynema_sgf_unit_test_exe_name} PRIVATE
  physics_test_utils.cpp

  test_simtime.cpp
  test_field.cpp
  test_field_ops.cpp
  test_field_fillpatch_ops.cpp
  test_physics.cpp
  test_auxiliary_fill.cpp
  )

add_subdirectory(vs)
