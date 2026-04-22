if (AMR_WIND_ENABLE_KYNEMA)
  target_sources(${kynema_sgf_lib_name} PRIVATE
    KynemaIface.cpp
    TurbineKynema.cpp
    )
endif()
