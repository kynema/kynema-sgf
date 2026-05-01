target_sources(${kynema_sgf_unit_test_exe_name} PRIVATE
  # testing utilities
  mms_test_utils.cpp

  # test cases
  test_mms_init.cpp
  test_mms_src.cpp
  test_mms_error.cpp
  )
