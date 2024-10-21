cmake_minimum_required (VERSION 3.20.1)

project(test-simulator-refinement)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)
  ## These test use dump diagnostics so require HighFive!
  phare_python3_exec(9  harris2d harris.py  ${CMAKE_CURRENT_BINARY_DIR})
  if(testMPI)
    phare_mpi_python3_exec(9 2  simple_2d_refinement test_2d_2_core.py  ${CMAKE_CURRENT_BINARY_DIR})
    phare_mpi_python3_exec(11 10 complex_2d_refinement test_2d_10_core.py  ${CMAKE_CURRENT_BINARY_DIR})
  endif(testMPI)

endif()
