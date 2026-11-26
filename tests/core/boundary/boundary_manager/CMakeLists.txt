cmake_minimum_required (VERSION 3.20.1)

project(test-boundary-manager)

function(_add_test test_name)

  add_executable(${test_name} ${test_name}.cpp)

  target_include_directories(${test_name} PRIVATE
    ${GTEST_INCLUDE_DIRS}
  )

  target_link_libraries(${test_name} PRIVATE
    phare_initializer
    phare_amr
    ${GTEST_LIBS}
)

  add_no_mpi_phare_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR})

endfunction()

_add_test(test_boundary_manager)
