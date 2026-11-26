cmake_minimum_required (VERSION 3.20.1)

project(test-conservation)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

# energy conservation vs regime: ~34 runs as parallel processes, each serial
# (tests/simulator/invariants.py has no MPI reduction), so never under mpirun.
# Level 11 keeps it out of pull-request CI; nightly label for the nightly pipeline
phare_mpi_python3_exec(11 1 conserv conserv.py ${CMAKE_CURRENT_BINARY_DIR})
phare_test_labels(conserv nightly)
