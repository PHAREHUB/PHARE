configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py @ONLY)

add_subdirectory(core)
add_subdirectory(initializer)
add_subdirectory(amr)
add_subdirectory(diagnostic)
add_subdirectory(simulator)
add_subdirectory(functional)
