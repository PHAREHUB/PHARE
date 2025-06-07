cmake_minimum_required (VERSION 3.20.1)

project(test-simulator-refinement)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)
  ## These test use dump diagnostics so require HighFive!
  if(testMPI)

    # disabled due to https://github.com/PHAREHUB/PHARE/issues/853
    # phare_mpi_python3_exec(11 10 complex_2d_refinement test_2d_10_core.py  ${CMAKE_CURRENT_BINARY_DIR})

  endif(testMPI)

  phare_python3_exec(9 simple_2d_refinement test_2d_2_core.py  ${CMAKE_CURRENT_BINARY_DIR})

endif()
