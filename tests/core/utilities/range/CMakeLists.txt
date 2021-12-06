cmake_minimum_required (VERSION 3.9)

project(test_range)

function(_add_range_test src_name)
  add_executable(${src_name} ${src_name}.cpp)
  target_include_directories(${src_name} PRIVATE ${GTEST_INCLUDE_DIRS} )
  target_link_libraries(${src_name} PRIVATE phare_core ${GTEST_LIBS})
  add_no_mpi_phare_test(${src_name} ${CMAKE_CURRENT_BINARY_DIR})
endfunction(_add_range_test)

_add_range_test(test_range)
_add_range_test(test_ranges)
_add_range_test(test_range_replacer)
