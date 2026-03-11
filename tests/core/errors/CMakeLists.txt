cmake_minimum_required (VERSION 3.20.1)

project(test-core)

function(_setup_core_test src_name)
  add_executable(${src_name} ${SOURCES_INC} ${src_name}.cpp)
  target_link_libraries(${src_name} PRIVATE phare_core cppdict ${GTEST_LIBS})
  add_no_mpi_phare_test(${src_name} ${CMAKE_CURRENT_BINARY_DIR})
endfunction(_setup_core_test)

_setup_core_test(test_fpe_errors)
target_compile_definitions(test_fpe_errors PRIVATE PHARE_CATCH_FPE)
