cmake_minimum_required (VERSION 3.3)

project(test-diagnostics)

set(SOURCES_INC
  ${PROJECT_SOURCE_DIR}/include.hpp
  ${PROJECT_SOURCE_DIR}/defs.hpp
  ${PROJECT_SOURCE_DIR}/integrator.hpp
  ${PROJECT_SOURCE_DIR}/samrai_lifecycle.hpp
  ${PROJECT_SOURCE_DIR}/tag_strat.hpp
  ${PROJECT_SOURCE_DIR}/hierarchy.hpp
)

set(SOURCES_CPP
  ${PROJECT_SOURCE_DIR}/test_main.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})

target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${gtest_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${gmock_SOURCE_DIR}/include>
  )

target_compile_options(${PROJECT_NAME} PRIVATE ${PHARE_WERROR_FLAGS})
target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_amr
  phare_solver
  gtest
  gmock)

add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME}
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

include(${PHARE_PROJECT_DIR}/sanitizer.cmake)
