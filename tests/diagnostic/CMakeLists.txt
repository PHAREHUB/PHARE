cmake_minimum_required (VERSION 3.9)

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

  if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    file(GLOB PYFILES "*.py")
    file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  message(STATUS "diagnostic working directory " ${PHARE_PROJECT_DIR})

endif()

