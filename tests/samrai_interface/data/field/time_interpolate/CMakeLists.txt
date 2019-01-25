cmake_minimum_required (VERSION 3.3)

project(test-field-data-time-interpolate)


set(SOURCES_INC
   )

set(SOURCES_CPP
  test_main.cpp
   )

add_executable(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})


target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${gtest_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${gmock_SOURCE_DIR}/include>
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_samrai_interface
  gtest
  gmock)


add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})


include(${PHARE_PROJECT_DIR}/sanitizer.cmake)
