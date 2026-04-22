if (AMR_WIND_ENABLE_OPENFAST)
  target_sources(${kynema_sgf_lib_name} PRIVATE
    FastIface.cpp
    TurbineFast.cpp
    )
endif()
