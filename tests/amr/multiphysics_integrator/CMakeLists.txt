cmake_minimum_required (VERSION 3.9)

project(test-multiphysics-integrator)

set(SOURCES_INC
   )

set(SOURCES_CPP
  test_main.cpp
   )

add_executable(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})


target_include_directories(${PROJECT_NAME} PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  ${GTEST_INCLUDE_DIRS}
  )

target_link_libraries(${PROJECT_NAME} PRIVATE
  phare_simulator
  pybind11::embed
  ${GTEST_LIBS})


add_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})



configure_file(${CMAKE_CURRENT_SOURCE_DIR}/job.py.in ${CMAKE_CURRENT_BINARY_DIR}/job.py @ONLY)
