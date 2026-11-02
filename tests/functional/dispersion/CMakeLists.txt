cmake_minimum_required (VERSION 3.20.1)

project(test-dispersion)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)
  ## These test use dump diagnostics so require HighFive!
  if(testMPI)
    phare_mpi_python3_exec(10 5 test-dispersion dispersion.py ${CMAKE_CURRENT_BINARY_DIR})
  else()
    phare_python3_exec(10 test-dispersion dispersion.py ${CMAKE_CURRENT_BINARY_DIR})
  endif()
  # 6 runs (3 interpolation orders x 2 resolutions), ~80 s in Release: runs on every PR
  phare_test_labels(test-dispersion physics)
endif()
