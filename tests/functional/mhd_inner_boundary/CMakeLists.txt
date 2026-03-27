cmake_minimum_required (VERSION 3.20.1)

project(test-mhd_inner_boundary)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)

  ## Visualization case: writes VTK output at t=0 to inspect inner boundary mesh data.
  phare_python3_exec(111, test-mhd_inner_boundary_sphere sphere_init.py ${CMAKE_CURRENT_BINARY_DIR})

endif()
