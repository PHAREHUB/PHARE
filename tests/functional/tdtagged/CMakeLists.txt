cmake_minimum_required (VERSION 3.20.1)

project(test-tdtagged)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)

  ## These test use dump diagnostics so require HighFive!
  # with MPI, 4 ranks rather than the default 2: level 0 then has patch boundaries
  # (x = 50, 150) that the sheets cross, which exercises tagging across patches
  if(testMPI)
    phare_mpi_python3_exec(11 4 test-td-tagged td1dtagged.py ${CMAKE_CURRENT_BINARY_DIR})
  else()
    phare_python3_exec(11 test-td-tagged td1dtagged.py ${CMAKE_CURRENT_BINARY_DIR})
  endif()

endif()
