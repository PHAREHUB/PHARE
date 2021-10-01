cmake_minimum_required (VERSION 3.9)

project(test-simulator-refinement)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)
  ## These test use dump diagnostics so require HighFive!
  if(testMPI)
    phare_mpi_python3_exec(9 2  simple_2d_refinement test_2d_2_core.py  ${CMAKE_CURRENT_BINARY_DIR})
    phare_mpi_python3_exec(11 10 complex_2d_refinement test_2d_10_core.py  ${CMAKE_CURRENT_BINARY_DIR})
  endif(testMPI)
endif()
