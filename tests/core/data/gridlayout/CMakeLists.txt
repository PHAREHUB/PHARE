cmake_minimum_required (VERSION 3.3)

project(test-GridLayout)

set(SOURCES_INC
  gridlayout_allocsize.h
  gridlayout_base_params.h
  gridlayout_cell_centered_coord.h
  gridlayout_deriv.h
  gridlayout_field_centered_coord.h
  gridlayout_indexing.h
  gridlayout_params.h
  gridlayout_test.h
  gridlayout_utilities.h
   )

set(SOURCES_CPP
  gridlayout_allocsize.cpp
  gridlayout_cell_centered_coord.cpp
  gridlayout_deriv.cpp
  gridlayout_field_centered_coord.cpp
  gridlayout_indexing.cpp
  test_main.cpp
   )
add_executable(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})



target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${gtest_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${gmock_SOURCE_DIR}/include>
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_core
  gtest
  gmock)

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(TARGET ${PROJECT_NAME}
  POST_BUILD
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gridlayout_generate_test_files.py
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

include(${PHARE_PROJECT_DIR}/sanitizer.cmake)
