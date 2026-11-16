cmake_minimum_required (VERSION 3.20.1)

project(test-fully-kinetic-waves)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)
  ## These test use dump diagnostics so require HighFive!
  # Langmuir, light wave and two-stream, one test per interpolation order (~2.5 min each, Debug)
  foreach(order 1 2 3)
    phare_python3_exec(10 test-fully-kinetic-waves-o${order} waves.py ${CMAKE_CURRENT_BINARY_DIR} ${order})
    phare_test_labels(test-fully-kinetic-waves-o${order} physics)
  endforeach()
endif()
