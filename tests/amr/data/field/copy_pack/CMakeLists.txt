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
  ${GTEST_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  )

target_link_libraries(${PROJECT_NAME} PUBLIC
  phare_amr
  ${GTEST_LIBS})

add_subdirectory(copy)
add_subdirectory(copy_overlap)
# add_subdirectory(stream_pack) # disabled. see https://github.com/PHAREHUB/PHARE/issues/60
