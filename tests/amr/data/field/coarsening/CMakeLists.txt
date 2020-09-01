cmake_minimum_required (VERSION 3.9)

project(test-field-data-coarsen)

configure_file(input_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/input_config.h)


set(SOURCES_INC
  test_basic_hierarchy.h
  test_linear_coarsen.h
  test_tag_strategy.h
  test_weighter_coarsen.h
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
  ${GTEST_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_amr
  ${GTEST_LIBS})


add_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(TARGET ${PROJECT_NAME}
                   POST_BUILD
                   COMMAND PYTHONPATH=${CMAKE_BINARY_DIR} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_coarsen_field.py ${CMAKE_CURRENT_BINARY_DIR})


