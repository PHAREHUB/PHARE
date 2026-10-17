cmake_minimum_required (VERSION 3.20.1)

project(test-diagnostics)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

set(SOURCES_INC
    ${PROJECT_SOURCE_DIR}/test_diagnostics.hpp
)

if(HighFive)
  # one binary per dim/interp, see SimulatorPermutations in tests/simulator/per_test.hpp
  #  to prevent OOM when compiling: all 2d permutations in one source file can exceed
  #  github actions runner memory (16GB) with -O3 -g
  function(_add_diagnostics_test dim interp)
    set(src_name test-diagnostics_${dim}d_interp${interp})

    add_executable(${src_name} ${SOURCES_INC} test-diagnostics.cpp)
    target_compile_definitions(${src_name} PRIVATE PHARE_TEST_DIM=${dim} PHARE_TEST_INTERP=${interp})

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

  foreach(dim 1 2 3)
    foreach(interp 1 2 3)
      _add_diagnostics_test(${dim} ${interp})
    endforeach()
  endforeach()

  message(STATUS "diagnostic working directory " ${PHARE_PROJECT_DIR})

endif()

