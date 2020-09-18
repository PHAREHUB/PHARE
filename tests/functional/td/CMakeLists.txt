cmake_minimum_required (VERSION 3.9)

project(test-td)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)

  ## These test use dump diagnostics so require HighFive!
  add_python3_test(td1d td1d.py ${CMAKE_CURRENT_BINARY_DIR})

endif()
