cmake_minimum_required (VERSION 3.20.1)

project(test-mhd-harris)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive AND testMPI)
  ## These test use dump diagnostics so require HighFive!
  # exec level 11
  # mpirun -n 10
  if(testMPI)
    phare_mpi_python3_exec(11 24 harris_2d harris_2d.py ${CMAKE_CURRENT_BINARY_DIR})
  endif(testMPI)
endif()


