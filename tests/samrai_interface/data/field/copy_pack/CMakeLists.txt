cmake_minimum_required (VERSION 3.3)

project(test-field-data)

set( SOURCES_INC
  field_data_test_param.h
  )
set ( SOURCES_CPP
    field_data_test_param.cpp
  )

add_library(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${gtest_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${gmock_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  )

target_link_libraries(${PROJECT_NAME} PUBLIC
  phare_samrai_interface
  gtest
  gmock)

add_subdirectory(copy)
add_subdirectory(copy_overlap)
add_subdirectory(stream_pack)

