cmake_minimum_required (VERSION 3.20.1)

project(test-field-data-stream-pack)


set(SOURCES_INC
  test_stream_pack_centered_ex.hpp
  test_stream_pack_centered_ey.hpp
   )

set(SOURCES_CPP
  test_main.cpp
  test_stream_pack_centered_ex.cpp
  test_stream_pack_centered_ey.cpp
   )

add_executable(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})

target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  test-field-data
  )


add_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})


