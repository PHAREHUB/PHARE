cmake_minimum_required (VERSION 3.20.1)

project(test-mhd_shock)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)

  ## These test use dump diagnostics so require HighFive!
	phare_mpi_python3_exec(21 10 test-mhd_shock mhd_shock.py ${CMAKE_CURRENT_BINARY_DIR})
endif()
