cmake_minimum_required (VERSION 3.3)

project(test-particle_selector)

set(SOURCES test_main.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${gtest_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${gmock_SOURCE_DIR}/include>
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_core
  gtest
  gmock)

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})

include(${PHARE_PROJECT_DIR}/sanitizer.cmake)
