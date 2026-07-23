cmake_minimum_required (VERSION 3.20.1)

project(test-simulator-amr-convergence)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)
  ## AMR space-time convergence studies (2D MHD); dump + read back h5 diagnostics.
  ## One ctest entry per test method so `ctest -j` runs the independent sweeps
  ## concurrently; inside an entry the sweep is sequential (order/drift need
  ## all points together).
  if(PHARE_PYTEST_SIMULATORS GREATER_EQUAL 2)
    phare_mpi_python3_exec(21 ${PHARE_MPI_PROCS} amr-convergence-alfven2d-spatial-o2      test_alfven2d_amr.py   ${CMAKE_CURRENT_BINARY_DIR} AlfvenConvergenceTest.test_spatial_convergence_1_2)
    phare_mpi_python3_exec(21 ${PHARE_MPI_PROCS} amr-convergence-alfven2d-spatial-legacy  test_alfven2d_amr.py   ${CMAKE_CURRENT_BINARY_DIR} AlfvenConvergenceTest.test_spatial_convergence_legacy)
    phare_mpi_python3_exec(21 ${PHARE_MPI_PROCS} amr-convergence-alfven2d-temporal        test_alfven2d_amr.py   ${CMAKE_CURRENT_BINARY_DIR} AlfvenConvergenceTest.test_temporal_sigma_sweep)
    phare_mpi_python3_exec(21 ${PHARE_MPI_PROCS} amr-convergence-whistler2d-spatial-o2     test_whistler2d_amr.py ${CMAKE_CURRENT_BINARY_DIR} WhistlerConvergenceTest.test_spatial_convergence_1_2)
    phare_mpi_python3_exec(21 ${PHARE_MPI_PROCS} amr-convergence-whistler2d-spatial-legacy test_whistler2d_amr.py ${CMAKE_CURRENT_BINARY_DIR} WhistlerConvergenceTest.test_spatial_convergence_legacy)
    phare_mpi_python3_exec(21 ${PHARE_MPI_PROCS} amr-convergence-whistler2d-temporal       test_whistler2d_amr.py ${CMAKE_CURRENT_BINARY_DIR} WhistlerConvergenceTest.test_temporal_sigma_sweep)
  endif()
endif()
