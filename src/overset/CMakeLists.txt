find_package(MPI REQUIRED COMPONENTS CXX)
find_package(TIOGA REQUIRED)
target_link_libraries_system(${kynema_sgf_lib_name} PUBLIC tioga)

target_sources(${kynema_sgf_lib_name}
  PRIVATE
  SGFTiogaIface.cpp
  TiogaInterface.cpp
  OversetOps.cpp
  overset_ops_routines.cpp
  )
