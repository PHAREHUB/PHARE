cmake_minimum_required (VERSION 3.20.1)

project(test-field-refine)

configure_file(input_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/input_config.h)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

set(SOURCES_INC
  test_tag_strategy.hpp
  ${CMAKE_CURRENT_BINARY_DIR}/input_config.h
   )

function(_setup_amr_field_refine_test src_name)

  add_executable(${src_name} ${SOURCES_INC} ${src_name}.cpp)

  target_include_directories(${src_name} PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    ${GTEST_INCLUDE_DIRS}
  )

  target_link_libraries(${src_name} PRIVATE phare_amr ${GTEST_LIBS})

  target_include_directories(${src_name} PRIVATE
    $<BUILD_INTERFACE:${SAMRAI_INCLUDE_DIRS}>)

  target_link_libraries(${src_name} PRIVATE ${SAMRAI_LIBRARIES})

endfunction(_setup_amr_field_refine_test)

function(_add_serial_amr_field_refine_test src_name)
  _setup_amr_field_refine_test(${src_name})
  add_no_mpi_phare_test(${src_name} ${CMAKE_CURRENT_BINARY_DIR})
endfunction(_add_serial_amr_field_refine_test)


_add_serial_amr_field_refine_test(test_default_field_refinement)

