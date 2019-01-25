cmake_minimum_required (VERSION 3.3)
project(test-field-refine)

configure_file(input_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/input_config.h)


set(SOURCES_INC
  test_basic_hierarchy.h
  test_tag_strategy.h
  ${CMAKE_CURRENT_BINARY_DIR}/input_config.h
   )

set(SOURCES_CPP
  test_basic_hierarchy.cpp
  test_main.cpp
  test_tag_strategy.cpp
   )

add_executable(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})

target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${gtest_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${gmock_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_samrai_interface
  gtest
  gmock)


target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${SAMRAI_INCLUDE_DIRS}>)

target_link_libraries(${PROJECT_NAME} PRIVATE ${SAMRAI_LIBRARIES})

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})

include(${PHARE_PROJECT_DIR}/sanitizer.cmake)
