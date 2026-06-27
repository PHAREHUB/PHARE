cmake_minimum_required (VERSION 3.20.1)

project(test-conservation)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

# These test use dump diagnostics so require HighFive!
if(HighFive)
  if(testMPI)
    phare_mpi_python3_exec(11 4 test-conservation conserv_2d.py ${CMAKE_CURRENT_BINARY_DIR})
  endif(testMPI)
endif()
