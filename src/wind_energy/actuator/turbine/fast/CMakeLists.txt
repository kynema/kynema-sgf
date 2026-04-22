if (AMR_WIND_ENABLE_OPENFAST)
  target_sources(${amr_wind_lib_name} PRIVATE
    FastIface.cpp
    TurbineFast.cpp
    )
endif()
