cmake_minimum_required (VERSION 3.3)

project(test-field-data-copy-with-overlap)


set(SOURCES_INC
  test_copy_overlap_centered_ex.h
  test_copy_overlap_centered_ey.h
   )

set(SOURCES_CPP
  test_copy_overlap_centered_ex.cpp
  test_copy_overlap_centered_ey.cpp
  test_main.cpp
   )

add_executable(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})

target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  test-field-data
  )


add_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})


