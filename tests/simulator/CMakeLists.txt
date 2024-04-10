cmake_minimum_required (VERSION 3.20.1)

project(test-simulator)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

phare_python3_exec(9 validation           test_validation.py      ${CMAKE_CURRENT_BINARY_DIR})
phare_python3_exec(9 data-wrangler        data_wrangler.py        ${CMAKE_CURRENT_BINARY_DIR})
phare_python3_exec(9 sim-refineParticlNbr refined_particle_nbr.py ${CMAKE_CURRENT_BINARY_DIR})
add_no_mpi_python3_test(periodicity test_init_periodicity.py ${CMAKE_CURRENT_BINARY_DIR})

if(HighFive)
  ## These test use dump diagnostics so require HighFive!
  phare_python3_exec(9       diagnostics test_diagnostics.py  ${CMAKE_CURRENT_BINARY_DIR}) # serial or n = 2
  phare_python3_exec(9       restarts    test_restarts.py     ${CMAKE_CURRENT_BINARY_DIR}) # serial or n = 2
  phare_python3_exec(9       run         test_run.py          ${CMAKE_CURRENT_BINARY_DIR}) # serial or n = 2
  # phare_python3_exec(9       tagging     test_tagging.py      ${CMAKE_CURRENT_BINARY_DIR}) # serial or n = 2

  if(testMPI)
    phare_mpi_python3_exec(9 3 diagnostics test_diagnostics.py  ${CMAKE_CURRENT_BINARY_DIR})
    phare_mpi_python3_exec(9 4 diagnostics test_diagnostics.py  ${CMAKE_CURRENT_BINARY_DIR})
  endif(testMPI)

  phare_python3_exec(11, test_diagnostic_timestamps test_diagnostic_timestamps.py ${CMAKE_CURRENT_BINARY_DIR})
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.py ${CMAKE_CURRENT_BINARY_DIR}/config.py @ONLY)

add_subdirectory(initialize)
add_subdirectory(refinement)
add_subdirectory(advance)

