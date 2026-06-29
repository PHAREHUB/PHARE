cmake_minimum_required (VERSION 3.20.1)

project(test-mhd_2d_vector_potential_init)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)

  ## This test uses dump diagnostics so requires HighFive!
  phare_python3_exec(21 test-mhd_2d_vector_potential_init test_vector_potential_init.py ${CMAKE_CURRENT_BINARY_DIR})
endif()
