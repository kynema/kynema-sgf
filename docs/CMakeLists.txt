# Setup paths
set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
set(DOXYGEN_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/doxygen)
set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/src)
set(DOXYGEN_OUTPUT_DIR ${SPHINX_SOURCE}/doxygen)
set(DOXYGEN_HTML_OUTPUT_DIR ${DOXYGEN_OUTPUT_DIR}/html)
file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})

# Doxygen target
find_package(Doxygen REQUIRED)

file(GLOB_RECURSE KYNEMA_SGF_PUBLIC_SOURCE ${DOXYGEN_INPUT_DIR}/*.H ${DOXYGEN_INPUT_DIR}/*.cpp)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_HTML_OUTPUT_DIR}/index.html)
set(DOXYGEN_TAG_FILE ${DOXYGEN_HTML_OUTPUT_DIR}/tagfile.xml)
set(DOXYGEN_AWESOME_STYLE_FILE ${PROJECT_SOURCE_DIR}/submods/doxygen-awesome-css/doxygen-awesome.css)
set(DOXYFILE_IN ${DOXYGEN_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${DOXYGEN_OUTPUT_DIR}/Doxyfile)
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                   DEPENDS ${KYNEMA_SGF_PUBLIC_SOURCE}
                   DEPENDS ${DOXYGEN_SOURCE_DIR}/main_page.md
                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                   MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN})

add_custom_target(doxygen DEPENDS ${DOXYGEN_INDEX_FILE}
                          COMMENT "Generating documentation with Doxygen")

# Doxysphinx target
find_package(Doxysphinx REQUIRED)

add_custom_target(doxysphinx
  COMMAND ${DOXYSPHINX_EXECUTABLE} build
  ${SPHINX_SOURCE} ${SPHINX_BUILD}/html ${DOXYGEN_HTML_OUTPUT_DIR}
  COMMENT "Converting doxygen documentation to sphinx with doxysphinx")
add_dependencies(doxysphinx doxygen)

# Sphinx target
find_package(Sphinx REQUIRED)

add_custom_target(sphinx-spelling
  COMMAND ${SPHINX_EXECUTABLE} -b spelling
  ${SPHINX_SOURCE} ${SPHINX_BUILD}/spelling
  COMMENT "Spell checking the documentation")
add_dependencies(sphinx-spelling doxysphinx)

add_custom_target(sphinx
  COMMAND ${SPHINX_EXECUTABLE} -M html
  ${SPHINX_SOURCE} ${SPHINX_BUILD}
  -c ${SPHINX_SOURCE} -W --keep-going -n
  COMMENT "Generating documentation with Sphinx")
add_dependencies(sphinx sphinx-spelling)

# All documentation
add_custom_target(docs ALL)
add_dependencies(docs sphinx)

