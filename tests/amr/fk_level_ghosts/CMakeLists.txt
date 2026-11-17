cmake_minimum_required (VERSION 3.20.1)

project(test-fk-level-ghosts)

# one binary per interpolation order, each with its own job (the order is compile-time)
foreach(order 1 2 3)
  set(binary ${PROJECT_NAME}-o${order})

  add_executable(${binary} test_main.cpp)
  target_compile_definitions(${binary} PRIVATE PHARE_TEST_INTERP_ORDER=${order})

  add_dependencies(${binary} cpp_etc)

  target_include_directories(${binary} PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    ${GTEST_INCLUDE_DIRS}
    )

  target_link_libraries(${binary} PRIVATE
    phare_simulator
    pybind11::embed
    ${GTEST_LIBS}
  )

  add_no_mpi_phare_test(${binary} ${CMAKE_CURRENT_BINARY_DIR})

  set(INTERP_ORDER ${order})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/job.py.in ${CMAKE_CURRENT_BINARY_DIR}/job_o${order}.py @ONLY)
endforeach()

# the level ghost band against every particle of the plasma, same jobs
foreach(order 1 2 3)
  set(binary test-fk-level-ghost-band-o${order})

  add_executable(${binary} test_band.cpp)
  target_compile_definitions(${binary} PRIVATE PHARE_TEST_INTERP_ORDER=${order})

  add_dependencies(${binary} cpp_etc)

  target_include_directories(${binary} PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    ${GTEST_INCLUDE_DIRS}
    )

  target_link_libraries(${binary} PRIVATE
    phare_simulator
    pybind11::embed
    ${GTEST_LIBS}
  )

  add_no_mpi_phare_test(${binary} ${CMAKE_CURRENT_BINARY_DIR})
endforeach()
