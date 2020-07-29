find_package(Doxygen)
if (NOT DOXYGEN_FOUND)
  message(FATAL_ERROR 
    "Doxygen is needed to build the documentation. Please install it correctly")
endif()

set(doc_inc_file ${PROJECT_SOURCE_DIR}/inc)
set(doc_src_file ${PROJECT_SOURCE_DIR}/src)
set(readme ${PROJECT_SOURCE_DIR}/README.md)
configure_file(Doxyfile
  ${PROJECT_BINARY_DIR}/Doxyfile  @ONLY IMMEDIATE)

add_custom_target(half_edge_doc ALL 
	COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
	SOURCES ${PROJECT_BINARY_DIR}/Doxyfile)
