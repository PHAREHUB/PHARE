cmake_minimum_required (VERSION 3.20.1)

project(test-simulator-init)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)
  ## These test use dump diagnostics so require HighFive!

  if(PHARE_PYTEST_SIMULATORS GREATER_EQUAL 1)
    phare_mpi_python3_exec(9 ${PHARE_MPI_PROCS} init-1d-fields     test_fields_init_1d.py     ${CMAKE_CURRENT_BINARY_DIR})
    phare_mpi_python3_exec(9 ${PHARE_MPI_PROCS} init-1d-particles  test_particles_init_1d.py  ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  if(PHARE_PYTEST_SIMULATORS GREATER_EQUAL 2)
    phare_mpi_python3_exec(9 ${PHARE_MPI_PROCS} init-2d-fields     test_fields_init_2d.py     ${CMAKE_CURRENT_BINARY_DIR})
    phare_mpi_python3_exec(9 ${PHARE_MPI_PROCS} init-2d-particles  test_particles_init_2d.py  ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  if(PHARE_PYTEST_SIMULATORS GREATER_EQUAL 3)
    phare_mpi_python3_exec(9 ${PHARE_MPI_PROCS} init-3d-fields     test_fields_init_3d.py     ${CMAKE_CURRENT_BINARY_DIR})
    phare_mpi_python3_exec(9 ${PHARE_MPI_PROCS} init-3d-particles  test_particles_init_3d.py  ${CMAKE_CURRENT_BINARY_DIR})
  endif()

  phare_mpi_python3_exec(9 ${PHARE_MPI_PROCS} init-densities density_check.py  ${CMAKE_CURRENT_BINARY_DIR})

endif()
