target_sources(${kynema_sgf_lib_name}
  PRIVATE

  ABL.cpp
  ABLStats.cpp
  ABLFieldInit.cpp
  ABLFieldInitFile.cpp
  ABLWallFunction.cpp
  MOData.cpp
  ABLMesoscaleForcing.cpp
  ABLMesoscaleInput.cpp
  ABLAnelastic.cpp
  )

add_subdirectory(actuator)
