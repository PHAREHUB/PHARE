cmake_minimum_required (VERSION 3.20.1)

project(test-tiles)

function(_add_tiles_test src_name)

  add_executable(${src_name} ${src_name}.cpp)

  target_include_directories(${src_name} PRIVATE
    ${GTEST_INCLUDE_DIRS}
  )

  target_link_libraries(${src_name} PRIVATE
    phare_core
    ${GTEST_LIBS})

  add_no_mpi_phare_test(${src_name} ${CMAKE_CURRENT_BINARY_DIR})

endfunction(_add_tiles_test)

_add_tiles_test(test_tile)
_add_tiles_test(test_tile_set_mapper)
