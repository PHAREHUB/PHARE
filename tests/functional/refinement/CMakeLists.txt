cmake_minimum_required (VERSION 3.20.1)

project(test-refinement)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive AND testMPI)
  ## divB preservation of the composite field-refinement kernels (order 2).
  ## These tests dump B diagnostics and read them back, so they require HighFive.
  # exec level 11
  phare_mpi_python3_exec(11 4 divb_refine_boxes   divb_refinement.py ${CMAKE_CURRENT_BINARY_DIR} boxes 2)
  phare_mpi_python3_exec(11 4 divb_refine_tagging divb_refinement.py ${CMAKE_CURRENT_BINARY_DIR} tagging 2)
endif()
