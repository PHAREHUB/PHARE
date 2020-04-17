cmake_minimum_required (VERSION 3.3)

project(test-diagnostics)

set(SOURCES_INC
    ${PROJECT_SOURCE_DIR}/test_diagnostics.h
)

set(SOURCES_CPP
  ${PROJECT_SOURCE_DIR}/test_diagnostics.cpp
)

if(HighFive)
  add_executable(${PROJECT_NAME} ${SOURCES_INC} ${SOURCES_CPP})

  target_include_directories(${PROJECT_NAME} PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    ${GTEST_INCLUDE_DIRS}
  )

  target_link_directories(${PROJECT_NAME} PRIVATE
    ${HDF5_LIBRARY_PATH}
  )

  target_link_libraries(${PROJECT_NAME} PRIVATE
    phare_simulator
    pybind11::embed
    ${GTEST_LIBS}
    ${HDF5_LIBRARIES})

  add_phare_test(${PROJECT_NAME} ${CMAKE_CURRENT_BINARY_DIR})

  message(STATUS "diagnostic working directory " ${PHARE_PROJECT_DIR})
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/job.py.in ${CMAKE_CURRENT_BINARY_DIR}/job.py @ONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py @ONLY)
