cmake_minimum_required (VERSION 3.20.1)

project(test-dispersion)

if(NOT ${PHARE_PROJECT_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  file(GLOB PYFILES "*.py")
  file(COPY ${PYFILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(HighFive)
  ## These test use dump diagnostics so require HighFive!
  # one test per interpolation order: the orders share nothing, so ctest runs them in
  # parallel (all three: 52 s -> ~20 s on 5 ranks in Release; ASan build was 26 min serial)
  foreach(order 1 2 3)
    if(testMPI)
      phare_mpi_python3_exec(10 5 test-dispersion-o${order} dispersion.py ${CMAKE_CURRENT_BINARY_DIR} ${order})
    else()
      phare_python3_exec(10 test-dispersion-o${order} dispersion.py ${CMAKE_CURRENT_BINARY_DIR} ${order})
    endif()
    phare_test_labels(test-dispersion-o${order} physics)
  endforeach()
endif()
