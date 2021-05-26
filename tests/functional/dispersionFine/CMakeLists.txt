cmake_minimum_required (VERSION 3.9)

project(test-dispersion-whistler-AIC)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)

  ## These test use dump diagnostics so require HighFive!
  phare_python3_exec(19, test-dispersion-whistler-AIC dispersion1d.py ${CMAKE_CURRENT_BINARY_DIR})
  phare_python3_exec(19, test-dispersion-whistler-AIC dispersion2d.py ${CMAKE_CURRENT_BINARY_DIR})
endif()
