
# target_link_libraries_system
#
# This function is similar to target_link_libraries but allows the includes
# determined from the library to be added as system includes to suppress
# warnings generated from those header files
#
# https://stackoverflow.com/questions/52135983/cmake-target-link-libraries-include-as-system-to-suppress-compiler-warnings/52136398#52136398
#
function(target_link_libraries_system target visibility)
  set(libs ${ARGN})
  foreach(lib ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${target} SYSTEM ${visibility} ${lib_include_dirs})
    target_link_libraries(${target} ${visibility} ${lib})
  endforeach(lib)
endfunction(target_link_libraries_system)

function(set_cuda_build_properties target)
  if (AMR_WIND_ENABLE_CUDA)
    get_target_property(_tgt_src ${target} SOURCES)
    list(FILTER _tgt_src INCLUDE REGEX "\\.cpp")
    set_source_files_properties(${_tgt_src} PROPERTIES LANGUAGE CUDA)
    set_target_properties(${target} PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
    if (AMR_WIND_ENABLE_CUDA_RDC)
      set_target_properties(${target} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    endif()
  endif()
endfunction(set_cuda_build_properties)

macro(init_amrex)
  if (${AMR_WIND_USE_INTERNAL_AMREX})
    set(AMREX_SUBMOD_LOCATION "${PROJECT_SOURCE_DIR}/submods/amrex")
    include(${PROJECT_SOURCE_DIR}/cmake/set_amrex_options.cmake)
    list(APPEND CMAKE_MODULE_PATH "${AMREX_SUBMOD_LOCATION}/Tools/CMake")
    add_subdirectory(${AMREX_SUBMOD_LOCATION})
    set(FCOMPARE_EXE ${CMAKE_BINARY_DIR}/submods/amrex/Tools/Plotfile/amrex_fcompare
      CACHE INTERNAL "Path to fcompare executable for regression tests")
  else()
    set(CMAKE_PREFIX_PATH ${AMREX_DIR} ${CMAKE_PREFIX_PATH})
    list(APPEND AMREX_COMPONENTS
      "3D" "PIC" "PARTICLES" "PDOUBLE" "DOUBLE" "LSOLVERS")
    if (AMR_WIND_ENABLE_MPI)
      list(APPEND AMREX_COMPONENTS "MPI")
    endif()
    if (AMR_WIND_ENABLE_OPENMP)
      list(APPEND AMREX_COMPONENTS "OMP")
    endif()
    if (AMR_WIND_ENABLE_CUDA)
      list(APPEND AMREX_COMPONENTS "CUDA")
    endif()
    if (AMR_WIND_ENABLE_SYCL)
      list(APPEND AMREX_COMPONENTS "SYCL")
    endif()
    if (AMR_WIND_ENABLE_ROCM)
      list(APPEND AMREX_COMPONENTS "HIP")
    endif()
    if (AMR_WIND_ENABLE_HYPRE)
      list(APPEND AMREX_COMPONENTS "HYPRE")
    endif()
    if (AMR_WIND_ENABLE_TINY_PROFILE)
      list(APPEND AMREX_COMPONENTS "TINY_PROFILE")
    endif()
    separate_arguments(AMREX_COMPONENTS)
    find_package(AMReX CONFIG REQUIRED
      COMPONENTS ${AMREX_COMPONENTS})
    message(STATUS "Found AMReX = ${AMReX_DIR}")
    set(FCOMPARE_EXE ${AMReX_DIR}/../../../bin/amrex_fcompare
      CACHE INTERNAL "Path to fcompare executable for regression tests")
  endif()
endmacro(init_amrex)

macro(init_amrex_hydro)
  if (${AMR_WIND_USE_INTERNAL_AMREX_HYDRO})
    set(AMREX_HYDRO_SUBMOD_LOCATION "${PROJECT_SOURCE_DIR}/submods/AMReX-Hydro")
    include(${PROJECT_SOURCE_DIR}/cmake/set_amrex_hydro_options.cmake)
    add_subdirectory(${AMREX_HYDRO_SUBMOD_LOCATION})
  else()
    set(CMAKE_PREFIX_PATH ${AMReX-Hydro_DIR} ${CMAKE_PREFIX_PATH})
    find_package(AMReX-Hydro CONFIG REQUIRED)
    message(STATUS "Found AMReX-Hydro = ${AMReX-Hydro_DIR}")
  endif()
endmacro(init_amrex_hydro)

macro(init_waves2amr)
  set(WAVES2AMR_SUBMOD_LOCATION "${PROJECT_SOURCE_DIR}/submods/Waves2AMR")
  include(${PROJECT_SOURCE_DIR}/cmake/set_waves2amr_options.cmake)
  add_subdirectory(${WAVES2AMR_SUBMOD_LOCATION})
endmacro(init_waves2amr)

macro(init_code_checks)
  if(AMR_WIND_ENABLE_CLANG_TIDY)
    set(CLANG_TIDY_EXEC_NAME "clang-tidy" CACHE STRING "Name of the clang-tidy executable")
    find_program(CLANG_TIDY_EXE NAMES "${CLANG_TIDY_EXEC_NAME}")
    if(CLANG_TIDY_EXE)
      message(STATUS "${CLANG_TIDY_EXEC_NAME} found: ${CLANG_TIDY_EXE}")
      set(CLANG_TIDY_EXE "${CLANG_TIDY_EXE};--config-file=${PROJECT_SOURCE_DIR}/.clang-tidy")
      find_program (CLANG_TIDY_CACHE_EXE NAMES "clang-tidy-cache")
      if(CLANG_TIDY_CACHE_EXE)
        message(STATUS "clang-tidy-cache found: ${CLANG_TIDY_CACHE_EXE}")
        set(CLANG_TIDY_EXE "${CLANG_TIDY_CACHE_EXE};${CLANG_TIDY_EXE}")
      endif()
    else()
      message(WARNING "${CLANG_TIDY_EXEC_NAME} not found.")
    endif()
  endif()

  if(AMR_WIND_ENABLE_CPPCHECK)
    find_program(CPPCHECK_EXE NAMES "cppcheck")
    if(CPPCHECK_EXE)
      message(STATUS "cppcheck found: ${CPPCHECK_EXE}")
      include(ProcessorCount)
      ProcessorCount(NP)
      if(NP EQUAL 0)
        set(NP 1)
      endif()
      file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/cppcheck)
      add_custom_target(cppcheck
          COMMAND ${CMAKE_COMMAND} -E echo "Running cppcheck on project using ${NP} cores..."
          COMMAND ${CMAKE_COMMAND} -E make_directory cppcheck
          # cppcheck ignores -isystem directories, so we change them to regular -I include directories (with no spaces either)
          COMMAND sed "s/isystem /I/g" ${CMAKE_BINARY_DIR}/compile_commands.json > cppcheck_compile_commands.json
          COMMAND ${CPPCHECK_EXE} --template=gcc --inline-suppr --suppress=unusedFunction --suppress=useStlAlgorithm --suppress=missingIncludeSystem --std=c++17 --language=c++ --enable=all --project=cppcheck_compile_commands.json -i ${PROJECT_SOURCE_DIR}/submods/amrex/Src -i ${PROJECT_SOURCE_DIR}/submods/AMReX-Hydro -i ${PROJECT_SOURCE_DIR}/submods/Waves2AMR -i ${PROJECT_SOURCE_DIR}/submods/googletest --output-file=cppcheck-full-report.txt -j ${NP}
          COMMENT "Run cppcheck on project compile_commands.json"
          BYPRODUCTS cppcheck-full-report.txt
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/cppcheck
          VERBATIM USES_TERMINAL
      )
      add_custom_target(cppcheck-ci
          # Filter out submodule source files after analysis
          COMMAND awk -v nlines=2 "/submods/ {for (i=0; i<nlines; i++) {getline}; next} 1" < cppcheck/cppcheck-full-report.txt > cppcheck/cppcheck-short-report.txt
          COMMAND cat cppcheck/cppcheck-short-report.txt | egrep "information:|error:|performance:|portability:|style:|warning:" | sort > cppcheck-ci-report.txt
          COMMAND printf "Warnings: " >> cppcheck-ci-report.txt
          COMMAND cat cppcheck-ci-report.txt | awk "END{print NR-1}" >> cppcheck-ci-report.txt
          COMMENT "Filter cppcheck results to only AMR-Wind files with results in cppcheck-ci-report.txt"
          DEPENDS cppcheck
          BYPRODUCTS cppcheck-ci-report.txt
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
          VERBATIM
      )

    else()
      message(WARNING "cppcheck not found.")
    endif()
  endif()
endmacro(init_code_checks)

macro(generate_version_info)
  include(GetGitRevisionDescription)
  get_git_head_revision(AMR_WIND_GITREFSPEC AMR_WIND_GIT_COMMIT_SHA)
  if (AMR_WIND_GIT_COMMIT_SHA)
    git_describe(AMR_WIND_VERSION_TAG "--tags" "--always")
    git_local_changes(AMR_WIND_REPO_DIRTY)
    option(AMR_WIND_HAVE_GIT_INFO "Git version for AMR-Wind" ON)
    if (${AMR_WIND_VERSION_TAG} MATCHES ".*-NOTFOUND")
      set(AMR_WIND_VERSION_TAG "v0.0.1")
    endif()
  endif()
  string(TIMESTAMP AMR_WIND_VERSION_TIMESTAMP "%Y-%m-%d %H:%M:%S (UTC)" UTC)
endmacro(generate_version_info)
