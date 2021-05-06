cmake_minimum_required (VERSION 3.9)

project(test-simulator-advance-1d)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)
  ## These test use dump diagnostics so require HighFive!

  add_python3_test(advance-1d-fields     fields.py     ${CMAKE_CURRENT_BINARY_DIR})
  add_python3_test(advance-1d-particles  particles.py  ${CMAKE_CURRENT_BINARY_DIR})


endif()
