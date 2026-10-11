cmake_minimum_required (VERSION 3.20.1)

project(test-dispersion)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)
  ## These test use dump diagnostics so require HighFive!
  if(testMPI)
    phare_mpi_python3_exec(11 5 test-dispersion dispersion.py ${CMAKE_CURRENT_BINARY_DIR})
  else()
    phare_python3_exec(11 test-dispersion dispersion.py ${CMAKE_CURRENT_BINARY_DIR})
  endif()
endif()
