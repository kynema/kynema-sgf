if (AMR_WIND_ENABLE_KYNEMA)
  target_sources(${amr_wind_lib_name} PRIVATE
    KynemaIface.cpp
    TurbineKynema.cpp
    )
endif()
