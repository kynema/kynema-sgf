target_sources(${amr_wind_lib_name}
  PRIVATE

  actuator_utils.cpp
  Actuator.cpp
  ActuatorContainer.cpp
  FLLC.cpp
  )

add_subdirectory(aero)
add_subdirectory(wing)
add_subdirectory(turbine)
add_subdirectory(disk)
