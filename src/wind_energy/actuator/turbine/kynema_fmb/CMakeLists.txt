if (KYNEMA_SGF_ENABLE_KYNEMA_FMB)
  target_sources(${kynema_sgf_lib_name} PRIVATE
    KynemaFMBIface.cpp
    TurbineKynemaFMB.cpp
    )
endif()
