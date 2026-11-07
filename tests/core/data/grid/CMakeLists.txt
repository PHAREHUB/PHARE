cmake_minimum_required (VERSION 3.20.1)

project(test-grid)


function(_add_grid_test src_name test_name)

  add_executable(${test_name} ${src_name}.cpp)

  target_include_directories(${test_name} PRIVATE
    ${GTEST_INCLUDE_DIRS}
  )

  target_link_directories(${test_name} PRIVATE
    ${HDF5_LIBRARY_PATH}
  )

  target_link_libraries(${test_name} PRIVATE
    phare_core
    phare_simulator
    pybind11::embed
    ${GTEST_LIBS})

  add_no_mpi_phare_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR})

endfunction(_add_grid_test)

_add_grid_test(test_grid test-grid)
_add_grid_test(test_grid_tile_set test_grid_tile_set)


