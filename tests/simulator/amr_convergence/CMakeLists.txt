cmake_minimum_required (VERSION 3.20.1)

project(test-simulator-amr-convergence)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)
  ## AMR space-time convergence studies (2D MHD); dump + read back h5 diagnostics.
  ## One ctest entry per file: the sweeps inside a file run sequentially, which is
  ## what the order and drift measurements need anyway (all points together).
  if(PHARE_PYTEST_SIMULATORS GREATER_EQUAL 2)
    phare_mpi_python3_exec(21 ${PHARE_MPI_PROCS} amr-convergence-alfven2d   test_alfven2d_amr.py   ${CMAKE_CURRENT_BINARY_DIR})
    phare_mpi_python3_exec(21 ${PHARE_MPI_PROCS} amr-convergence-whistler2d test_whistler2d_amr.py ${CMAKE_CURRENT_BINARY_DIR})
  endif()
endif()
