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
  # one binary per simulator permutation to prevent OOM when compiling/linking: several
  #  simulators in one binary can exceed github actions runner memory (16GB) with -O3 -g,
  #  or with coverage instrumentation
  function(_add_diagnostics_test dim interp refined_parts)
    set(src_name test-diagnostics_${dim}d_interp${interp}_refined${refined_parts})

    add_executable(${src_name} ${SOURCES_INC} test-diagnostics.cpp)
    target_compile_definitions(${src_name} PRIVATE
      PHARE_TEST_DIM=${dim} PHARE_TEST_INTERP=${interp} PHARE_TEST_REFINED_PARTS=${refined_parts})

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

  set(_refined_parts_1_1 2 3)
  set(_refined_parts_1_2 2 3 4)
  set(_refined_parts_1_3 2 3 4 5)

  set(_refined_parts_2_1 4 5 8 9)
  set(_refined_parts_2_2 4 5 8 9 16)
  set(_refined_parts_2_3 4 5 8 25)

  set(_refined_parts_3_1 6)
  set(_refined_parts_3_2 6)
  set(_refined_parts_3_3 6)

  foreach(dim 1 2 3)
    foreach(interp 1 2 3)
      foreach(refined_parts ${_refined_parts_${dim}_${interp}})
        _add_diagnostics_test(${dim} ${interp} ${refined_parts})
      endforeach()
    endforeach()
  endforeach()

  message(STATUS "diagnostic working directory " ${PHARE_PROJECT_DIR})

endif()

