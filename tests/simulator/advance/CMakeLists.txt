cmake_minimum_required (VERSION 3.9)

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

  endif()

  add_no_mpi_python3_test(multithreading_nd test_multithreading.py ${CMAKE_CURRENT_BINARY_DIR})

endif()
