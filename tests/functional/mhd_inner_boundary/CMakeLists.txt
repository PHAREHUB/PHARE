cmake_minimum_required (VERSION 3.20.1)

project(test-mhd-inner-boundary)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)

  ## These tests use dump diagnostics so require HighFive!
  phare_python3_exec(101, test-mhd-inner-boundary-sphere sphere_periodic.py ${CMAKE_CURRENT_BINARY_DIR})
endif()
