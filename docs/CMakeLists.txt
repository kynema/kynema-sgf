# Setup paths
set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx)
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)

# Sphinx target
find_package(Sphinx REQUIRED)

add_custom_target(sphinx-spelling
  COMMAND ${SPHINX_EXECUTABLE} -b spelling
  ${SPHINX_SOURCE} ${SPHINX_BUILD}/spelling
  COMMENT "Spell checking the documentation")

add_custom_target(sphinx
  COMMAND ${SPHINX_EXECUTABLE} -M html
  ${SPHINX_SOURCE} ${SPHINX_BUILD}
  -c ${SPHINX_SOURCE} -W --keep-going -n
  COMMENT "Generating documentation with Sphinx")
add_dependencies(sphinx sphinx-spelling)

# All documentation
add_custom_target(docs ALL)
add_dependencies(docs sphinx)

