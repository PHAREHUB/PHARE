cmake_minimum_required (VERSION 3.3)

project(test-resource)

set(SOURCES_INC
  basic_hierarchy.h
  cell_placeholder.h
  resource_test_1d.h
   )

set(SOURCES_CPP
  basic_hierarchy.cpp
  resource_test_1d.cpp
  test_main.cpp
   )

add_executable(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})

target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${gtest_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${gmock_SOURCE_DIR}/include>
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_amr
  gtest
  gmock)


target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${SAMRAI_INCLUDE_DIRS}>)

target_link_libraries(${PROJECT_NAME} PRIVATE ${SAMRAI_LIBRARIES})

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR})

include(${PHARE_PROJECT_DIR}/sanitizer.cmake)
