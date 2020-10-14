cmake_minimum_required (VERSION 3.9)

project(test-particles-data-refine)

configure_file(input_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/input_config.h)

set(SOURCES_INC
  test_particledata_refine_basic_hierarchy.h
  test_particle_data_refine_tag_strategy.h
  ${CMAKE_CURRENT_BINARY_DIR}/input_config.h
   )

set(SOURCES_CPP
  test_particle_data_refine.cpp
  test_split.cpp
   )

add_executable(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})


target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  ${GTEST_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_amr
  ${GTEST_LIBS})


add_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})


