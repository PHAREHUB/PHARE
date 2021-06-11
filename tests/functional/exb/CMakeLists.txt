cmake_minimum_required (VERSION 3.9)

project(test-exb2d)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)

  ## These test use dump diagnostics so require HighFive!
  phare_python3_exec(9, test-exb2d exb.py ${CMAKE_CURRENT_BINARY_DIR})
endif()
