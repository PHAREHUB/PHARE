cmake_minimum_required (VERSION 3.20.1)

project(test-diagnostics)

set(SOURCES_INC
    ${PROJECT_SOURCE_DIR}/test_diagnostics.hpp
)

if(HighFive)
  function(_add_diagnostics_test src_name)

    add_executable(${src_name} ${SOURCES_INC} ${src_name}.cpp)

    target_include_directories(${src_name} PRIVATE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
      ${GTEST_INCLUDE_DIRS}
    )

    target_link_directories(${src_name} PRIVATE
      ${HDF5_LIBRARY_PATH}
    )

    target_link_libraries(${src_name} PRIVATE
      phare_simulator
      pybind11::embed
      ${GTEST_LIBS}
      ${HDF5_LIBRARIES})

    add_dependencies(${src_name} cpp_etc)

    add_phare_test(${src_name} ${CMAKE_CURRENT_BINARY_DIR})

  endfunction(_add_diagnostics_test)

  _add_diagnostics_test(test-diagnostics_1d)
  _add_diagnostics_test(test-diagnostics_2d)
  _add_diagnostics_test(test-diagnostics_3d)

  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/job_1d.py.in ${CMAKE_CURRENT_BINARY_DIR}/job_1d.py @ONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/job_2d.py.in ${CMAKE_CURRENT_BINARY_DIR}/job_2d.py @ONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/job_3d.py.in ${CMAKE_CURRENT_BINARY_DIR}/job_3d.py @ONLY)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py @ONLY)

  message(STATUS "diagnostic working directory " ${PHARE_PROJECT_DIR})

endif()

