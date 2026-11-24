target_sources(${kynema_sgf_lib_name}
  PRIVATE

  actuator_utils.cpp
  Actuator.cpp
  ActuatorContainer.cpp
  FLLC.cpp
  )

add_subdirectory(aero)
add_subdirectory(motion)
add_subdirectory(wing)
add_subdirectory(sector)
add_subdirectory(drone)
add_subdirectory(turbine)
add_subdirectory(disk)
