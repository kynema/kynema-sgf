target_sources(${kynema_sgf_unit_test_exe_name}
  PRIVATE
  utest_main.cpp
  test_config.cpp
  )

add_subdirectory(ks_test_utils)
add_subdirectory(core)
add_subdirectory(wind_energy)
add_subdirectory(offshore_wind)
add_subdirectory(utilities)
add_subdirectory(equation_systems)
add_subdirectory(turbulence)
add_subdirectory(fvm)
add_subdirectory(multiphase)
add_subdirectory(ocean_waves)
add_subdirectory(projection)
add_subdirectory(boundary_conditions)
add_subdirectory(eos_models)

if(KYNEMA_SGF_ENABLE_MASA)
  add_subdirectory(mms)
endif()

target_compile_options(${kynema_sgf_unit_test_exe_name} PRIVATE
  $<$<COMPILE_LANGUAGE:CXX>:${KYNEMA_SGF_CXX_FLAGS}>)
target_include_directories(${kynema_sgf_unit_test_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

#target_link_libraries(${utest_exe_name} PRIVATE GTest::GTest GTest::Main)
target_link_libraries(${kynema_sgf_unit_test_exe_name} PRIVATE gtest)
target_include_directories(${kynema_sgf_unit_test_exe_name} SYSTEM PRIVATE
  ${PROJECT_SOURCE_DIR}/submods/googletest/googletest/include)

target_link_libraries(${kynema_sgf_unit_test_exe_name} PRIVATE ${kynema_sgf_lib_name} AMReX-Hydro::amrex_hydro_api)
if (KYNEMA_SGF_ENABLE_W2A)
  target_link_libraries(${kynema_sgf_unit_test_exe_name} PRIVATE Waves2AMR::Waves2AMR)
endif()

install(TARGETS ${kynema_sgf_unit_test_exe_name})
