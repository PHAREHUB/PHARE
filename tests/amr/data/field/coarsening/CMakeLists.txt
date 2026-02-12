cmake_minimum_required (VERSION 3.20.1)

project(test-field-data-coarsen)


set(SOURCES_INC
  test_weighter_coarsen.hpp
  )

set(SOURCES_CPP test_coarsening_main.cpp)

add_executable(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})


target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  ${GTEST_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_amr
  phare_initializer
  ${GTEST_LIBS})


add_no_mpi_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})
