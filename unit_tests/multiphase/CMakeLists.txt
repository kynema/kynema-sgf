target_sources(
  ${kynema_sgf_unit_test_exe_name} PRIVATE
  test_channel_builder.cpp
  test_vof_plic.cpp
  test_vof_cons.cpp
  test_vof_tools.cpp
  test_momflux.cpp
  test_vof_BCs.cpp
  test_mflux_schemes.cpp
  test_reference_fields.cpp
  test_vof_overset_ops.cpp
  )
