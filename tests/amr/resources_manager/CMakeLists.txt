cmake_minimum_required (VERSION 3.3)

project(test-resource)

configure_file(input_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/input_config.h)

set(SOURCES_INC
  basic_hierarchy.h
  ${CMAKE_CURRENT_BINARY_DIR}/input_config.h
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
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  )

target_compile_options(${PROJECT_NAME} PRIVATE ${PHARE_WERROR_FLAGS})
target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_amr
  phare_initializer
  gtest
  gmock)


add_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})

include(${PHARE_PROJECT_DIR}/sanitizer.cmake)
