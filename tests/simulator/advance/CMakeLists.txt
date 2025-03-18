cmake_minimum_required (VERSION 3.20.1)

project(test-simulator-advance)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)
  ## These test use dump diagnostics so require HighFive!

  phare_mpi_python3_exec(9 ${PHARE_MPI_PROCS} advance-1d-fields     test_fields_advance_1d.py     ${CMAKE_CURRENT_BINARY_DIR})
  phare_mpi_python3_exec(9 ${PHARE_MPI_PROCS} advance-1d-particles  test_particles_advance_1d.py  ${CMAKE_CURRENT_BINARY_DIR})

  if(NOT lowResourceTests)

    phare_mpi_python3_exec(9 ${PHARE_MPI_PROCS} advance-2d-fields     test_fields_advance_2d.py     ${CMAKE_CURRENT_BINARY_DIR})
    phare_mpi_python3_exec(9 ${PHARE_MPI_PROCS} advance-2d-particles  test_particles_advance_2d.py  ${CMAKE_CURRENT_BINARY_DIR})

    if(highResourceTests) # off by default as it's quite intensive
        phare_mpi_python3_exec(9 ${PHARE_MPI_PROCS} advance-3d-fields     test_fields_advance_3d.py     ${CMAKE_CURRENT_BINARY_DIR})
        phare_mpi_python3_exec(9 ${PHARE_MPI_PROCS} advance-3d-particles  test_particles_advance_3d.py  ${CMAKE_CURRENT_BINARY_DIR})
    endif(highResourceTests)

  endif()

endif()
