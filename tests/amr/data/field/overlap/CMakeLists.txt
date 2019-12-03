cmake_minimum_required (VERSION 3.3)

project(test-field-overlap)

# configure_file(input_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/input_config.h)

set(SOURCES_INC
  # ${CMAKE_CURRENT_BINARY_DIR}/input_config.h
   )

set(SOURCES_CPP
  test_main.cpp
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


target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${SAMRAI_INCLUDE_DIRS}>)

target_link_libraries(${PROJECT_NAME} PRIVATE ${SAMRAI_LIBRARIES})

add_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})


