target_sources(${kynema_sgf_unit_test_exe_name} PRIVATE
  test_actuator_sampling.cpp
  test_actuator_flat_plate.cpp
  test_actuator_fixed_wing.cpp
  test_airfoil.cpp
  test_actuator_free_functions.cpp
  test_disk_uniform_ct.cpp
  test_FLLC.cpp
  test_actuator_joukowsky_disk.cpp
  test_disk_functions.cpp
  )

if (KYNEMA_SGF_ENABLE_OPENFAST)
  target_sources(${kynema_sgf_unit_test_exe_name} PRIVATE
    test_fast_iface.cpp
    test_turbine_fast.cpp
    )
endif()
